clc; clear; close all;

% 1. Setup and Loading
fs_target = 16000;
file_path = 'data/clean_speech.wav';

if ~exist('outputs/plots', 'dir'), mkdir('outputs/plots'); end
[x, fs] = audioread(file_path);

if fs ~= fs_target
    x = resample(x, fs_target, fs);
    fs = fs_target;
end

t = (0:length(x)-1)/fs;

% 2. Figure 1: Time Domain Waveform
figure('Name', 'Reference Waveform', 'Color', 'w', 'Position', [100, 100, 800, 400]);
plot(t, x, 'b');
xlabel('Time (s)'); ylabel('Amplitude');
title('REFERENCE: Clean Speech Waveform');
grid on; axis tight; ylim([-1.1 1.1]);
saveas(gcf, 'outputs/plots/REF_WAVEFORM.png');

% 3. Figure 2: Magnitude Spectrum (Dynamic Centering)
L = length(x);
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f_axis = fs*(0:(floor(L/2)))/L;

mag_db = 20*log10(P1 + eps); % Convert to dB scale

% Dynamic Centering Logic:
% Find the peak value and set limits relative to it
peak_val = max(mag_db); 
y_min = peak_val - 60;  % Show 60 dB range below the peak
y_max = peak_val + 10;  % Leave 10 dB headroom above the peak

figure('Name', 'Reference Spectrum', 'Color', 'w', 'Position', [150, 150, 800, 400]);
plot(f_axis/1000, mag_db, 'b');
xlabel('Frequency (kHz)'); ylabel('Mag (dB)');
title('REFERENCE: Clean Speech Spectrum');
grid on; xlim([0 8]); 
ylim([y_min y_max]); % Apply dynamic limits to center the signal
saveas(gcf, 'outputs/plots/REF_SPECTRUM.png');

% 4. Figure 3: Spectrogram Analysis (Dynamic Contrast)
figure('Name', 'Reference Spectrogram', 'Color', 'w', 'Position', [200, 200, 800, 400]);
[s, f_spec, t_spec] = spectrogram(x, hamming(512), 256, 512, fs, 'yaxis');
spec_db = 20*log10(abs(s) + eps);

% Find peak in spectrogram for dynamic color contrast
max_spec = max(spec_db(:));

imagesc(t_spec, f_spec/1000, spec_db);
axis xy; colormap jet; 
caxis([max_spec-60 max_spec]); % Center color scale based on peak power
xlabel('Time (s)'); ylabel('Freq (kHz)');
title('REFERENCE: Clean Speech Spectrogram');
colorbar;
saveas(gcf, 'outputs/plots/REF_SPECTROGRAM.png');

fprintf('Reference analysis complete. Dynamic plots saved to outputs/plots/ folder.\n');