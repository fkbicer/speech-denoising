clc; clear; close all;

% 1. SETUP & LOADING
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_impulsive_5dB.wav'; 

% Ensure output directories exist
if ~exist('outputs/wavs', 'dir'), mkdir('outputs/wavs'); end
if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Input files not found! Please check the data/ directory.');
end

[x_ref_orig, fs_r] = audioread(file_clean);
[x_noiz_orig, fs_n] = audioread(file_noisy);

% Resample to project standard
x_ref = resample(x_ref_orig, fs_target, fs_r);
x_noiz = resample(x_noiz_orig, fs_target, fs_n);

% Align signal lengths
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== IMPULSIVE NOISE DENOISING ===\n');

% 2. IMPULSE DETECTION (Statistical Analysis)
% We calculate the moving standard deviation to identify sudden spikes
moving_std = movstd(abs(x_n), round(fs_target * 0.01)); % 10ms window
threshold = 3.5 * median(moving_std); % Sensitivity threshold

% Create a mask for samples exceeding the local statistical average
impulse_mask = abs(x_n) > (threshold + mean(abs(x_n)));

% 3. ADAPTIVE FILTERING
% Use a 5-point median filter as a replacement source
y_med = medfilt1(x_n, 5); 

% Adaptive Replacement: Only change samples where an impulse was detected
% This prevents the median filter from blurring the entire speech signal
y_proc = x_n;
y_proc(impulse_mask) = y_med(impulse_mask);

% 4. POST-PROCESSING (Smoothing)
% Apply a 6th order Low-pass filter at 7kHz to remove residual switching noise
[b, a] = butter(6, 7000/(fs_target/2), 'low');
y_final = filtfilt(b, a, y_proc);

% Normalization to prevent clipping
y_final = real(y_final);
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% 5. METRICS & VISUALIZATION
% Calculate SNR using the provided project formula
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Results (Impulsive Noise):\n');
fprintf('  Input SNR : %.2f dB\n', snr_pre);
fprintf('  Output SNR: %.2f dB\n', snr_post);
fprintf('  SNR Gain  : +%.2f dB\n', snr_post - snr_pre);

audiowrite('outputs/wavs/denoised_impulsive.wav', y_final, fs_target);
generate_report_plots(x_c, x_n, y_final, fs_target, 'Impulsive Noise');
saveas(gcf, 'outputs/plots/report_impulsive.png');