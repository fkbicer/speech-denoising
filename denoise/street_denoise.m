clc; clear; close all;

% 1. SETUP & LOADING (Task 1)
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_street_5dB.wav'; 

% Create output directories if they do not exist
if ~exist('outputs/wavs', 'dir'), mkdir('outputs/wavs'); end
if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Input files missing! Please verify the data/ folder.');
end

[x_ref_orig, fs_r] = audioread(file_clean);
[x_noiz_orig, fs_n] = audioread(file_noisy);

% Resample signals to 16 kHz project standard
x_ref = resample(x_ref_orig, fs_target, fs_r);
x_noiz = resample(x_noiz_orig, fs_target, fs_n);

% Synchronize signal lengths
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== STREET NOISE DENOISING (Task 3 & 4) ===\n');

% 2. PRE-PROCESSING: RUMBLE SUPPRESSION
% Street noise contains heavy low-frequency components from engines and tires.
% A 4th order high-pass filter at 200Hz targets this initial rumble.
[b_lc, a_lc] = butter(4, 200/(fs_target/2), 'high');
x_pre = filtfilt(b_lc, a_lc, x_n);

% 3. STRATEGY: SMOOTHED SPECTRAL ATTENUATION (Task 3)
win_len = 512; hop_len = 256; nfft = 512;
win = hamming(win_len);

[S, f_stft, t_stft] = stft(x_pre, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
Mag = abs(S);
Ph  = angle(S);

% Time-wise spectral smoothing to handle non-stationary traffic noise
Mag_smoothed = movmean(Mag, 3, 2); 

% Frequency-dependent attenuation
% We apply heavier suppression on frequencies where traffic noise dominates
freqs = linspace(0, fs_target/2, size(Mag, 1))';
attenuation_mask = ones(size(freqs));
attenuation_mask(freqs < 400) = 0.4;   % Target residual engine drone
attenuation_mask(freqs > 4500) = 0.7;  % Target tire and wind hiss

Mag_denoised = Mag_smoothed .* attenuation_mask;

% 4. ENERGY-BASED GATING
% Soft gate to suppress background street noise during speech pauses
gate_threshold = 0.05 * max(Mag_denoised(:));
gate_mask = Mag_denoised > gate_threshold;
Mag_denoised = Mag_denoised .* gate_mask;

% 5. RECONSTRUCTION (ISTFT)
S_new = Mag_denoised .* exp(1j * Ph);
y_raw = istft(S_new, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
y_raw = real(y_raw);

% Critical Fix: Length Alignment
if length(y_raw) >= L
    y_denoised = y_raw(1:L);
else
    y_denoised = [y_raw; zeros(L - length(y_raw), 1)];
end

% 6. POST-PROCESSING: CLARITY BOOST
% Add a subtle high-frequency lift to restore speech intelligibility
[b_air, a_air] = butter(2, 3500/(fs_target/2), 'high');
air_band = filtfilt(b_air, a_air, y_denoised);
y_final = y_denoised + 0.15 * air_band;

% Final Normalization
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% 7. METRICS & REPORTING
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Results (Street Noise):\n');
fprintf('  Input SNR : %.2f dB\n', snr_pre);
fprintf('  Output SNR: %.2f dB\n', snr_post);
fprintf('  SNR Gain  : +%.2f dB\n', snr_post - snr_pre);

audiowrite('outputs/wavs/denoised_street.wav', y_final, fs_target);
generate_report_plots(x_c, x_n, y_final, fs_target, 'Street Noise');
saveas(gcf, 'outputs/plots/report_street.png');