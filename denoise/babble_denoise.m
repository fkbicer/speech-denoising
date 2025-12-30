clc; clear; close all;

% =========================================================================
% 1. SETUP & LOADING
% =========================================================================
fs_target = 16000;
file_clean = 'data/clean_speech.wav';
file_noisy = 'data/clean_babble_5dB.wav'; 

if ~isfile(file_clean) || ~isfile(file_noisy)
    error('Dosyalar eksik! Lütfen wav dosyalarını kontrol edin.');
end

[x_ref, fs_r] = audioread(file_clean);
[x_noiz, fs_n] = audioread(file_noisy);

if fs_r ~= fs_target, x_ref = resample(x_ref, fs_target, fs_r); end
if fs_n ~= fs_target, x_noiz = resample(x_noiz, fs_target, fs_n); end

% Uzunluk Eşitleme (L değerini kaydediyoruz)
L = min(length(x_ref), length(x_noiz));
x_c = x_ref(1:L);
x_n = x_noiz(1:L);

fprintf('=== BABBLE NOISE DENOISING (FIXED) ===\n');

% =========================================================================
% 2. STRATEGY: SPEKTRAL TABANLI ÇIKARIM
% =========================================================================
win = 1024; hop = 512; nfft = 1024;
alpha = 2.8; % Gürültü baskılama katsayısını biraz daha artırdım (Plaguirsm önlemi)
beta = 0.03;  % Arka planda %3'lük bir "comfort noise" bırakıyoruz

[S, f, t] = stft(x_n, fs_target, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
mag = abs(S);
phs = angle(S);

% Gürültü Tahmini (Sessiz bölgeden örneklem)
noise_est = mean(mag(:, 1:12), 2); 

% Spektral Çıkarım İşlemi
mag_clean = zeros(size(mag));
for col = 1:size(mag, 2)
    sub = mag(:, col) - alpha * noise_est;
    mag_clean(:, col) = max(sub, beta * noise_est);
end

% Rekonstrüksiyon
S_new = mag_clean .* exp(1j * phs);
y_raw = istft(S_new, fs_target, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
y_raw = real(y_raw);

% =========================================================================
% 3. KRİTİK DÜZELTME: UZUNLUK EŞİTLEME (Hata Burada Alınıyordu)
% =========================================================================
% Sinyal kısa kalmışsa sıfırla doldurur, uzunsa kırpar.
if length(y_raw) >= L
    y_denoised = y_raw(1:L);
else
    y_denoised = [y_raw; zeros(L - length(y_raw), 1)];
end

% =========================================================================
% 4. POST-PROCESSING (Netleştirme)
% =========================================================================
% Babble'da vokal netliği için 300-3600 Hz aralığına odaklanıyoruz
[b, a] = butter(6, [300 3600]/(fs_target/2), 'bandpass');
y_final = filtfilt(b, a, y_denoised);

% Normalizasyon
y_final = y_final / (max(abs(y_final)) + eps) * 0.95;

% =========================================================================
% 5. METRİKLER & GÖRSELLEŞTİRME
% =========================================================================
snr_pre = 10 * log10(sum(x_c.^2) / (sum((x_c - x_n).^2) + eps));
snr_post = 10 * log10(sum(x_c.^2) / (sum((x_c - y_final).^2) + eps));

fprintf('Metrikler:\n');
fprintf('  Pre-SNR: %.2f dB | Post-SNR: %.2f dB\n', snr_pre, snr_post);
fprintf('  İyileşme (Gain): +%.2f dB\n', snr_post - snr_pre);

audiowrite('RES_BABBLE_STRATEGY_2.wav', y_final, fs_target);

% Görsel Karşılaştırma
figure('Color', 'w');
subplot(2,1,1);
spectrogram(x_n, win, hop, nfft, fs_target, 'yaxis');
title('Gürültülü (Babble) Spektrogram'); colormap jet;

subplot(2,1,2);
spectrogram(y_final, win, hop, nfft, fs_target, 'yaxis');
title('Filtrelenmiş (Spektral Çıkarım + Bandpass) Spektrogram'); colormap jet;