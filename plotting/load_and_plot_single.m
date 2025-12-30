clc;
clear;
close all;

%% Load clean speech
[x_clean, fs] = audioread('data/clean_speech.wav');

fprintf('Sampling frequency: %d Hz\n', fs);
fprintf('Signal length: %.2f seconds\n', length(x_clean)/fs);

%% Time-domain plot
t = (0:length(x_clean)-1)/fs;

figure;
plot(t, x_clean);
xlabel('Time (s)');
ylabel('Amplitude');
title('Clean Speech - Time Domain');
