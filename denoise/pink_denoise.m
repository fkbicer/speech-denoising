clc; clear; close all;

% 1. SETUP & LOADING (Task 1)
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_pink_5dB.wav'; 

% Ensure output directories exist
if ~exist('outputs/wavs', 'dir'), mkdir('outputs/wavs'); end
if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Input files not found! Please check the data/ directory.');
end

[x_ref_orig, fs_r] = audioread(file_clean);
[x_noiz_orig, fs_n] = audioread(file_noisy);

% Resample signals to project standard (16 kHz)
x_ref = resample(x_ref_orig, fs_target, fs_r);
x_noiz = resample(x_noiz_orig, fs_target, fs_n);

% Synchronize signal lengths
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== PINK NOISE DENOISING (Task 3 & 4) ===\n');

% 2. PRE-PROCESSING: PRE-EMPHASIS
% Pink noise has a 1/f characteristic with heavy low-frequency energy.
% Pre-emphasis flattens the spectrum by boosting high frequencies.
pre_emph_coeff = [1, -0.97];
x_pre = filter(pre_emph_coeff, 1, x_n);

% 3. STRATEGY: WIENER FILTERING IN STFT DOMAIN (Task 3)
win_len = 512; hop_len = 256; nfft = 1024;
win = hamming(win_len, 'periodic');

[S, f, t] = stft(x_pre, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
mag_noisy = abs(S);
phase_noisy = angle(S);

% Noise Estimation from the initial silent period
noise_est = mean(mag_noisy(:, 1:12), 2); 

% Wiener Gain Calculation
mag_sq = mag_noisy.^2;
noise_sq = noise_est.^2;
gain_factor = 1.6; % Aggressiveness factor for Pink noise

% Compute the Wiener gain function: G = P_sig / (P_sig + P_noise)
gain = max(0.02, (mag_sq - gain_factor * noise_sq) ./ (mag_sq + eps));
gain = movmedian(gain, 3, 2); % Spectral smoothing to reduce artifacts

% Apply gain and transform back to time domain
mag_denoised = mag_noisy .* gain;
S_clean = mag_denoised .* exp(1j * phase_noisy);
x_res = istft(S_clean, fs_target, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);

% 4. POST-PROCESSING & RECONSTRUCTION
y_raw = real(x_res);

% De-emphasis: Restore the natural tonal balance of the speech
y_deemp = filter(1, pre_emph_coeff, y_raw);

% Critical Fix: Length Alignment
if length(y_deemp) >= L
    y_final = y_deemp(1:L);
else
    y_final = [y_deemp; zeros(L - length(y_deemp), 1)];
end

% 150Hz High-pass filter to remove residual low-end rumble
[b_hp, a_hp] = butter(4, 150/(fs_target/2), 'high');
y_final = filtfilt(b_hp, a_hp, y_final);

% Normalize the output to prevent clipping
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% 5. METRICS & REPORTING (Task 2, 5, 6)
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Results (Pink Noise):\n');
fprintf('  Input SNR : %.2f dB\n', snr_pre);
fprintf('  Output SNR: %.2f dB\n', snr_post);
fprintf('  SNR Gain  : +%.2f dB\n', snr_post - snr_pre);

% Save the cleaned audio file
audiowrite('outputs/wavs/denoised_pink.wav', y_final, fs_target);

% Generate comparative report plots
generate_report_plots(x_c, x_n, y_final, fs_target, 'Pink Noise');

% Save the report figure
saveas(gcf, 'outputs/plots/report_pink.png');

fprintf('\n>>> Process Complete: Output files saved to outputs/wavs/ and outputs/plots/\n');