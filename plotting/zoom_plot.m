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

win_sec = 0.05;               
start_sec = 1.00;               

for k = 1:numel(files)
    [x, fs] = audioread(files{k});
    Nwin = round(win_sec * fs);
    i0 = round(start_sec * fs) + 1;
    i1 = min(i0 + Nwin - 1, length(x));

    t = (0:(i1-i0))/fs;

    figure;
    plot(t, x(i0:i1));
    grid on;
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Zoomed Time Domain (', num2str(win_sec*1000), ' ms) - ' names{k}]);

    ymax = max(abs(x(i0:i1)));
    if ymax > 0
        ylim([-1.1*ymax 1.1*ymax]);
    end
end
