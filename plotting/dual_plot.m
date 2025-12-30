clc; clear; close all;

[xc, fs] = audioread('data/clean_speech.wav');
[xn, ~ ] = audioread('data/clean_white_5dB.wav');

start_sec = 1.0;
win_sec   = 0.02;
i0 = round(start_sec*fs)+1;
i1 = i0 + round(win_sec*fs)-1;

t = (0:(i1-i0))/fs;

figure;
plot(t, xc(i0:i1)); hold on;
plot(t, xn(i0:i1));
grid on;
legend('Clean','Noisy');
title('Clean vs Noisy (same 20ms segment)');
xlabel('Time (s)'); ylabel('Amplitude');
