clc; clear; close all;


% 1. SETUP & LOADING
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_babble_5dB.wav'; 

% Ensure the required folders exist for organized output
if ~exist('outputs/wavs', 'dir'), mkdir('outputs/wavs'); end
if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Files missing! Please check the data/ folder.');
end

[x_ref_orig, fs_r] = audioread(file_clean);
[x_noiz_orig, fs_n] = audioread(file_noisy);

% Resampling to project standard (16 kHz)
x_ref = resample(x_ref_orig, fs_target, fs_r);
x_noiz = resample(x_noiz_orig, fs_target, fs_n);

% Length Alignment (Matching clean and noisy signals)
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== BABBLE NOISE DENOISING (Task 3 & 4) ===\n');

% 2. STRATEGY: SPECTRAL SUBTRACTION
% Babble noise is non-stationary and overlaps with speech frequencies.
% Spectral subtraction in the STFT domain is used to suppress background talkers.
win = 1024; hop = 512; nfft = 1024;
alpha = 2.8; % Over-subtraction factor (Aggressiveness of noise removal)
beta = 0.03;  % Spectral floor (Prevents "musical noise" by leaving a small floor)

[S, f, t] = stft(x_n, fs_target, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
mag = abs(S);
phs = angle(S);

% Noise Estimation (Using the first few frames as silence/noise reference)
noise_est = mean(mag(:, 1:12), 2); 

% Processing Spectral Magnitude
mag_clean = zeros(size(mag));
for col = 1:size(mag, 2)
    % Subtract estimated noise profile from the signal
    sub = mag(:, col) - alpha * noise_est;
    
    % Apply spectral floor to avoid zeroing out bins (minimizes artifacts)
    mag_clean(:, col) = max(sub, beta * noise_est);
end

% Reconstruction (ISTFT)
S_new = mag_clean .* exp(1j * phs);
y_raw = istft(S_new, fs_target, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
y_raw = real(y_raw);

% 3. CRITICAL FIX: LENGTH MATCHING
% Pad with zeros if short, or trim if longer than original L
if length(y_raw) >= L
    y_denoised = y_raw(1:L);
else
    y_denoised = [y_raw; zeros(L - length(y_raw), 1)];
end

% 4. POST-PROCESSING (Speech Enhancement)
% Apply a Butterworth Bandpass filter to isolate the vocal range (300-3600 Hz)
[b, a] = butter(6, [300 3600]/(fs_target/2), 'bandpass');
y_final = filtfilt(b, a, y_denoised);

% Normalization (Prevents clipping in the output file)
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% 5. METRICS & REPORTING
% Quantitative performance metrics using SNR formula provided in project
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Results (Babble Noise):\n');
fprintf('  Input SNR : %.2f dB\n', snr_pre);
fprintf('  Output SNR: %.2f dB\n', snr_post);
fprintf('  SNR Gain  : +%.2f dB\n', snr_post - snr_pre);

audiowrite('outputs/wavs/denoised_babble.wav', y_final, fs_target);
generate_report_plots(x_c, x_n, y_final, fs_target, 'Babble Noise');
saveas(gcf, 'outputs/plots/report_babble.png');