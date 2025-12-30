clc; clear; close all;

% 1. SETUP & LOADING (Task 1)
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_white_5dB.wav'; 

% Create output directories for better organization
if ~exist('outputs/wavs', 'dir'), mkdir('outputs/wavs'); end
if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Files missing! Check the data/ folder.');
end

[x_ref_orig, fs_r] = audioread(file_clean);
[x_noiz_orig, fs_n] = audioread(file_noisy);

% Resample to project standard 16 kHz
x_ref = resample(x_ref_orig, fs_target, fs_r);
x_noiz = resample(x_noiz_orig, fs_target, fs_n);

% Align signal lengths
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== WHITE NOISE DENOISING (Task 3 & 4) ===\n');

% 2. PRE-FILTERING: HISS SUPPRESSION
% White noise energy is constant across all frequencies.
% A 6th order Butterworth LPF at 5.5kHz removes the harsh high-frequency hiss.
[b_lp, a_lp] = butter(6, 5500/(fs_target/2), 'low');
x_pre = filtfilt(b_lp, a_lp, x_n);

% 3. STRATEGY: WIENER FILTERING IN STFT DOMAIN (Task 3)
win_len = 512; hop_len = 256; nfft = 512;
win = hamming(win_len);

[S, f_stft, t_stft] = stft(x_pre, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
Mag = abs(S);
Ph  = angle(S);

% Stationary Noise Estimation
% Since white noise is stationary, we estimate noise power across all frames
noise_power = median(Mag.^2, 2); 

% Compute Wiener Gain: G = P_sig / (P_sig + P_noise)
Gain = Mag.^2 ./ (Mag.^2 + noise_power + eps);
Mag_wiener = Mag .* Gain;

% 4. SPECTRAL GATING
% Soft thresholding to suppress residual noise floor during speech pauses
gate_threshold = 0.06 * max(Mag_wiener(:));
gate_mask = Mag_wiener > gate_threshold;
Mag_denoised = Mag_wiener .* gate_mask;

% 5. RECONSTRUCTION (ISTFT)
S_new = Mag_denoised .* exp(1j * Ph);
y_raw = istft(S_new, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
y_raw = real(y_raw);

% Critical Fix: Length Alignment
if length(y_raw) >= L
    y_final = y_raw(1:L);
else
    y_final = [y_raw; zeros(L - length(y_raw), 1)];
end

% 6. POST-PROCESSING: CLARITY RECOVERY
% High-frequency "Air" boost to restore natural speech brilliance after LPF
[b_air, a_air] = butter(2, 4000/(fs_target/2), 'high');
air_signal = filtfilt(b_air, a_air, y_final);
y_final = y_final + 0.15 * air_signal;

% Final Normalization
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% 7. METRICS & REPORTING (Task 2, 5, 6)
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Results (White Noise):\n');
fprintf('  Input SNR : %.2f dB\n', snr_pre);
fprintf('  Output SNR: %.2f dB\n', snr_post);
fprintf('  SNR Gain  : +%.2f dB\n', snr_post - snr_pre);

audiowrite('outputs/wavs/denoised_white.wav', y_final, fs_target);
generate_report_plots(x_c, x_n, y_final, fs_target, 'White Noise');
saveas(gcf, 'outputs/plots/report_white.png');