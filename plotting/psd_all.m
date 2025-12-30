clc; clear; close all;

files = {
    'data/clean_speech.wav'
    'data/clean_white_5dB.wav'
    'data/clean_pink_5dB.wav'
    'data/clean_impulsive_5dB.wav'
    'data/clean_babble_5dB.wav'
    'data/clean_street_5dB.wav'
};

names = {
    'Clean'
    'White (5 dB)'
    'Pink (5 dB)'
    'Impulsive (5 dB)'
    'Babble (5 dB)'
    'Street (5 dB)'
};
[xc, fs] = audioread(files{1});
Nfft = 65536;  % power of 2, fixed
f = (0:(Nfft/2)) * (fs/Nfft);

Xc = fft(xc, Nfft);
MagC = abs(Xc(1:Nfft/2+1));
MagC_dB = 20*log10(MagC + 1e-12);

for k = 2:numel(files)
    [x, fs2] = audioread(files{k});
    if fs2 ~= fs
        error('Sampling rate mismatch!');
    end

    X = fft(x, Nfft);
    Mag = abs(X(1:Nfft/2+1));
    Mag_dB = 20*log10(Mag + 1e-12);

    figure;
    plot(f, MagC_dB); hold on;
    plot(f, Mag_dB);
    grid on;
    xlim([0 8000]); 
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(['FFT Magnitude: Clean vs ' names{k}]);
    legend('Clean', names{k}, 'Location','best');
end
