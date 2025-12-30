clc; clear;

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
    'White'
    'Pink'
    'Impulsive'
    'Babble'
    'Street'
};

thr = 0.999;

for k = 1:numel(files)
    [x, fs] = audioread(files{k});
    clip_ratio = mean(abs(x) >= thr) * 100;
    fprintf('%-10s | fs=%d | clip=%.2f %%\n', names{k}, fs, clip_ratio);
end
