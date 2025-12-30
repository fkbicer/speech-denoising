%% street_denoise.m  (Optimized Street Noise Strategy)
clc; clear; close all;

%% Paths
in_clean = 'data/clean_speech.wav';
in_noisy = 'data/clean_street_5dB.wav';

scriptDir = fileparts(mfilename('fullpath'));
outdir = fullfile(scriptDir,'outputs');
if exist(outdir,'file')==2
    error('"outputs" exists as a FILE.');
end
if ~exist(outdir,'dir'), mkdir(outdir); end

out_wav = fullfile(outdir,'denoised_street.wav');

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
fprintf('STREET | Input SNR: %.2f dB\n', snr_in);

%% ===================== DENOISING PIPELINE =====================

%% 1) Adaptive rumble suppression (low-cut)
[b_lc, a_lc] = butter(4, 200/(fs/2), 'high');
x_lc = filtfilt(b_lc, a_lc, x);

%% 2) STFT
win  = hamming(512);
hop  = 256;
nfft = 512;

[S,~,~] = stft(x_lc, fs, ...
    'Window', win, ...
    'OverlapLength', hop, ...
    'FFTLength', nfft);

Mag = abs(S);
Ph  = angle(S);

%% 3) Spectral smoothing (street ≠ stationary)
Mag_s = movmean(Mag, 3, 2);   % time-wise smoothing

%% 4) Frequency-dependent attenuation
freq = linspace(0, fs/2, size(Mag_s,1))';
atten = ones(size(freq));

atten(freq < 300) = 0.3;     % suppress residual rumble
atten(freq > 5000) = 0.8;    % mild HF attenuation

Mag_d = Mag_s .* atten;

%% 5) Soft energy-based gate
thresh = 0.05 * max(Mag_d(:));
mask = Mag_d > thresh;
Mag_d = Mag_d .* mask;

%% 6) Reconstruction
S_new = Mag_d .* exp(1j*Ph);
y_tmp = istft(S_new, fs, ...
    'Window', win, ...
    'OverlapLength', hop, ...
    'FFTLength', nfft);

y_tmp = real(y_tmp);

%% 7) Length safety
if length(y_tmp) > N
    y_tmp = y_tmp(1:N);
elseif length(y_tmp) < N
    y_tmp = [y_tmp; zeros(N-length(y_tmp),1)];
end

%% 8) Gentle clarity boost
[b_air,a_air] = butter(2, 3500/(fs/2), 'high');
air = filtfilt(b_air, a_air, y_tmp);
y = y_tmp + 0.12*air;

%% 9) Normalize
mx = max(abs(y))+1e-12;
y = y/mx*0.95;

%% ===================== METRICS =====================
snr_out = snr_db(s,y);
fprintf('STREET | Output SNR: %.2f dB | Gain: %+0.2f dB\n', ...
    snr_out, snr_out-snr_in);

%% ===================== EXPORT =====================
audiowrite(out_wav, y, fs);
disp(['Saved: ' out_wav]);

%% ===================== PLOTS =====================
plot_analysis_pair(s, x, fs, 'STREET - Before');
plot_analysis_pair(s, y, fs, 'STREET - After');;

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
