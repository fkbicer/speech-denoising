%% white_denoise.m  (Optimized White Noise Strategy)
clc; clear; close all;

%% Paths
in_clean = 'data/clean_speech.wav';
in_noisy = 'data/clean_white_5dB.wav';

scriptDir = fileparts(mfilename('fullpath'));
outdir = fullfile(scriptDir, 'outputs');
if exist(outdir,'file')==2
    error('"outputs" exists as a FILE.');
end
if ~exist(outdir,'dir'), mkdir(outdir); end

out_wav = fullfile(outdir,'denoised_white.wav');

%% Load
[s, fs]  = audioread(in_clean); s = s(:,1);
[x, fs2] = audioread(in_noisy); x = x(:,1);
assert(fs==fs2,'Sampling rate mismatch');

N = min(length(s),length(x));
s = s(1:N); 
x = x(1:N);

%% SNR
snr_db = @(c,y) 10*log10(sum(c.^2)/(sum((y-c).^2)+1e-12));
snr_in = snr_db(s,x);

fprintf('WHITE | Input SNR: %.2f dB\n', snr_in);

%% ===================== DENOISING PIPELINE =====================

%% 1) Pre-filter: LPF to suppress hiss
[b_lp,a_lp] = butter(6, 5500/(fs/2), 'low');
x_lpf = filtfilt(b_lp, a_lp, x);

%% 2) STFT Analysis
win  = hamming(512);
hop  = 256;
nfft = 512;

[S,~,~] = stft(x_lpf, fs, ...
    'Window', win, ...
    'OverlapLength', hop, ...
    'FFTLength', nfft);

Mag = abs(S);
Ph  = angle(S);

%% 3) Wiener filtering in magnitude domain
% ---- Noise power estimate (time-wise minimum) ----
noise_power = median(Mag.^2, 2);   % stationary white noise assumption

% ---- Wiener gain ----
Gain = Mag.^2 ./ (Mag.^2 + noise_power);

% ---- Apply Wiener filtering ----
Mag_w = Mag .* Gain;

%% 4) Soft spectral gating
thresh = 0.06 * max(Mag_w(:));
mask = Mag_w > thresh;
Mag_d = Mag_w .* mask;

%% 5) Reconstruct (phase preserved)
S_new = Mag_d .* exp(1j*Ph);
y_tmp = istft(S_new, fs, ...
    'Window', win, ...
    'OverlapLength', hop, ...
    'FFTLength', nfft);

y_tmp = real(y_tmp);

%% 6) Air boost (clarity recovery)
[b_air,a_air] = butter(2, 4000/(fs/2), 'high');
air = filtfilt(b_air, a_air, y_tmp);
y = y_tmp + 0.15*air;

%% 7) Length + normalize
if length(y) > N
    y = y(1:N);
elseif length(y) < N
    y = [y; zeros(N - length(y), 1)];
end
mx = max(abs(y))+1e-12;
y = y/mx*0.95;

%% ===================== METRICS =====================
snr_out = snr_db(s,y);
fprintf('WHITE | Output SNR: %.2f dB | Gain: %+0.2f dB\n', ...
    snr_out, snr_out-snr_in);

%% ===================== EXPORT =====================
audiowrite(out_wav, y, fs);
disp(['Saved: ' out_wav]);

%% ===================== PLOTS =====================
plot_analysis_pair(s, x, fs, 'WHITE - Before');
plot_analysis_pair(s, y, fs, 'WHITE - After');

%% ==================================================
%% ===================== Local function =====================
function plot_analysis_pair(clean, sig, fs, tag)
    % Time zoom
    start_sec = 1.0;
    win_sec   = 0.02; % 20 ms

    i0 = round(start_sec*fs) + 1;
    i0 = max(1, min(i0, length(clean)));
    i1 = min(i0 + round(win_sec*fs) - 1, length(clean));
    t = (0:(i1-i0))/fs;

    figure('Name',['TimeZoom ' tag]);
    plot(t, clean(i0:i1), 'LineWidth', 1.2); hold on;
    plot(t, sig(i0:i1),   'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title(['Clean vs Signal (' num2str(win_sec*1000) ' ms) - ' tag]);
    legend('Clean','Signal','Location','best');

    % FFT magnitude
    Nfft = 65536;
    f = (0:Nfft/2) * (fs/Nfft);

    C = fft(clean, Nfft);
    S = fft(sig,   Nfft);

    MagC = 20*log10(abs(C(1:Nfft/2+1)) + 1e-12);
    MagS = 20*log10(abs(S(1:Nfft/2+1)) + 1e-12);

    figure('Name',['FFT ' tag]);
    plot(f, MagC); hold on;
    plot(f, MagS);
    grid on; xlim([0 fs/2]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(['FFT Magnitude - ' tag]);
    legend('Clean','Signal','Location','best');

    % Spectrogram
    figure('Name',['Spec ' tag]);
    subplot(2,1,1);
    spectrogram(clean, 512, 384, 1024, fs, 'yaxis');
    title(['Spectrogram - Clean (' tag ')']);

    subplot(2,1,2);
    spectrogram(sig, 512, 384, 1024, fs, 'yaxis');
    title(['Spectrogram - Signal (' tag ')']);
end
