function denoise_pink()
    clc; clear; close all;

    % 1. AYARLAR VE DOSYA OKUMA
    file_noisy = 'data/clean_pink_5dB.wav'; 
    file_clean_ref = 'data/clean_speech.wav'; 
    fs_target = 16000;

    if ~isfile(file_noisy), error('Gürültülü dosya bulunamadı!'); end
    [x_noisy, fs] = audioread(file_noisy);
    if fs ~= fs_target, x_noisy = resample(x_noisy, fs_target, fs); fs = fs_target; end

    % 2. ÖN İŞLEME (Pre-emphasis)
    pre_emph = [1, -0.97];
    x_pre = filter(pre_emph, 1, x_noisy);

    % 3. STFT PARAMETRELERİ
    win_len = 512;
    hop_len = 256;
    nfft = 1024;
    win = hamming(win_len, 'periodic');

    % Spektrogram hesapla
    [S, f, t] = stft(x_pre, fs, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);
    mag_noisy = abs(S);
    phase_noisy = angle(S);

    % 4. GÜRÜLTÜ TAHMİNİ
    noise_frames = round(0.3 * fs / hop_len);
    noise_mu = mean(mag_noisy(:, 1:noise_frames), 2); 

    % 5. GELİŞTİRİLMİŞ WIENER FİLTRELEME
    mag_sq = mag_noisy.^2;
    noise_sq = noise_mu.^2;
    
    % Gain (Kazanç) hesaplama - Gürültü baskılama faktörünü 1.5 olarak bıraktık
    gain = max(0.01, (mag_sq - 1.5 * noise_sq) ./ mag_sq); 
    gain = movmedian(gain, 3, 2); 

    % Uygulama
    mag_denoised = mag_noisy .* gain;

    % 6. REKONSTRÜKSİYON (ISTFT)
    S_clean = mag_denoised .* exp(1j * phase_noisy);
    x_res = istft(S_clean, fs, 'Window', win, 'OverlapLength', win_len - hop_len, 'FFTLength', nfft);

    % 7. SON İŞLEME
    % ÖNEMLİ: ISTFT sonrası oluşabilecek hayali kısımları temizle
    y_denoised = real(x_res); 
    
    % De-emphasis (Ön işlemi geri al)
    y_denoised = filter(1, pre_emph, y_denoised);
    
    % Rumble cut (Düşük frekans temizliği)
    [b, a] = butter(4, 150/(fs/2), 'high');
    y_denoised = filtfilt(b, a, y_denoised);

    % Normalizasyon
    y_denoised = y_denoised / (max(abs(y_denoised)) + eps) * 0.95;

    % 8. ANALİZ VE KAYIT
    if isfile(file_clean_ref)
        [x_c, fs_c] = audioread(file_clean_ref);
        if fs_c ~= fs, x_c = resample(x_c, fs, fs_c); end
        L = min([length(x_c), length(x_noisy), length(y_denoised)]);
        
        snr_pre = 10 * log10(sum(x_c(1:L).^2) / (sum((x_c(1:L) - x_noisy(1:L)).^2) + eps));
        snr_post = 10 * log10(sum(x_c(1:L).^2) / (sum((x_c(1:L) - y_denoised(1:L)).^2) + eps));
        
        fprintf('\nPink Noise Denoising Tamamlandı\n');
        fprintf('-------------------------------\n');
        fprintf('Giriş SNR: %.2f dB\n', snr_pre);
        fprintf('Çıkış SNR: %.2f dB\n', snr_post);
        fprintf('Kazanç:    +%.2f dB\n', snr_post - snr_pre);
    end

    % Artık hata vermeyecektir
    % Klasör kontrolü ve kayıt
    if ~exist('outputs', 'dir'), mkdir('outputs'); end
    output_name = 'outputs/PINK_CLEANED_FINAL.wav';
    
    audiowrite(output_name, y_denoised, fs);
    fprintf('\nSonuç buraya kaydedildi: %s\n', output_name);
    
    % Görselleştirme
    figure('Color', 'w', 'Name', 'Pink Noise Denoising Result');
    subplot(2,1,1);
    imagesc(t, f/1000, 20*log10(mag_noisy + eps)); axis xy; title('Noisy Spectrogram (Pink)');
    colorbar; caxis([-80 0]); ylabel('Freq (kHz)');
    
    subplot(2,1,2);
    imagesc(t, f/1000, 20*log10(mag_denoised + eps)); axis xy; title('Cleaned Spectrogram');
    colorbar; caxis([-80 0]); ylabel('Freq (kHz)'); xlabel('Time (s)');
end