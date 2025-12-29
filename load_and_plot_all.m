clc;
clear;
close all;

%% File list (dataset)
files = {
    'data/clean_speech.wav'
    'data/clean_white_5dB.wav'
    'data/clean_pink_5dB.wav'
    'data/clean_impulsive_5dB.wav'
    'data/clean_babble_5dB.wav'
    'data/clean_street_5dB.wav'
};

titles = {
    'Clean Speech'
    'White Noise (5 dB)'
    'Pink Noise (5 dB)'
    'Impulsive Noise (5 dB)'
    'Babble Noise (5 dB)'
    'Street Noise (5 dB)'
};

%% Loop over all signals
for k = 1:length(files)

    % Load signal
    [x, fs] = audioread(files{k});
    t = (0:length(x)-1)/fs;

    % Plot
    figure;
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Time Domain - ' titles{k}]);

    fprintf('%s loaded | Duration: %.2f sec\n', ...
        titles{k}, length(x)/fs);
end
