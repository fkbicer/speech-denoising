function denoise_impulsive()
    clc; clear; close all;

    % 1. SETUP & LOADING
    file_noisy = 'data/clean_impulsive_5dB.wav'; 
    file_clean_ref = 'data/clean_speech.wav'; 
    output_folder = 'outputs';
    fs_target = 16000;

    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    if ~isfile(file_noisy), error('Gürültülü dosya bulunamadı!'); end

    [x_n, fs] = audioread(file_noisy);
    if fs ~= fs_target, x_n = resample(x_n, fs_target, fs); fs = fs_target; end

    % 2. IMPULSE DETECTION (İstatistiksel Tespit)
    % Sinyalin mutlak değerinin yerel ortalamadan sapmasına bakıyoruz
    moving_std = movstd(abs(x_n), round(fs*0.01)); % 10ms pencere
    threshold = 3.5 * median(moving_std); % Eşik değer (Hassasiyet için 3.5 idealdir)
    
    % Darbeli noktaların maskesini oluştur
    impulse_mask = abs(x_n) > (threshold + mean(abs(x_n)));

    % 3. ADAPTIVE FILTERING (Sadece Tespit Edilen Yerleri Onar)
    x_proc = x_n;
    
    % Medyan filtre penceresini biraz genişletiyoruz (5 örnek)
    y_med = medfilt1(x_n, 5); 
    
    % Sadece maskelenmiş (gürültülü) kısımları medyan filtrelenmiş haliyle değiştir
    x_proc(impulse_mask) = y_med(impulse_mask);

    % 4. POST-PROCESSING (Pürüzsüzleştirme)
    % Filtreleme sonrası oluşan ufak çıtlamaları yok etmek için çok hafif LPF
    [b, a] = butter(6, 7000/(fs/2), 'low');
    y_denoised = filtfilt(b, a, x_proc);

    % ISTFT hatalarını önlemek için reel kısım garantisi ve normalizasyon
    y_denoised = real(y_denoised);
    y_denoised = y_denoised / (max(abs(y_denoised)) + eps) * 0.95;

    % 5. METRICS CALCULATION
    if isfile(file_clean_ref)
        [x_c, fs_c] = audioread(file_clean_ref);
        if fs_c ~= fs, x_c = resample(x_c, fs, fs_c); end
        L = min(length(x_c), length(y_denoised));
        
        % Standart SNR
        snr_pre = 10 * log10(sum(x_c(1:L).^2) / (sum((x_c(1:L) - x_n(1:L)).^2) + eps));
        snr_post = 10 * log10(sum(x_c(1:L).^2) / (sum((x_c(1:L) - y_denoised(1:L)).^2) + eps));
        
        % Segmental SNR (Konuşma kalitesi için daha tutarlıdır)
        seg_snr = calc_seg_snr(x_c(1:L), y_denoised(1:L), fs);

        fprintf('\nImpulsive Noise Denoising Tamamlandı\n');
        fprintf('------------------------------------\n');
        fprintf('Giriş SNR:      %.2f dB\n', snr_pre);
        fprintf('Çıkış SNR:     %.2f dB\n', snr_post);
        fprintf('Segmental SNR: %.2f dB\n', seg_snr);
        fprintf('Toplam Kazanç: +%.2f dB\n', snr_post - snr_pre);
    end

    % 6. SAVING & PLOTTING
    out_path = fullfile(output_folder, 'IMPULSIVE_CLEANED.wav');
    audiowrite(out_path, y_denoised, fs);
    
    figure('Color', 'w', 'Position', [100 100 900 600]);
    t = (0:length(y_denoised)-1)/fs;
    
    subplot(2,1,1);
    plot(t, x_n, 'r'); hold on;
    plot(t(impulse_mask), x_n(impulse_mask), 'k.'); % Tespit edilen darbeler
    title('Gürültülü Sinyal (Siyah noktalar tespit edilen darbelerdir)');
    ylabel('Amplitude'); axis tight; grid on;

    subplot(2,1,2);
    plot(t, y_denoised, 'g');
    title('Temizlenmiş Sinyal (Adaptive Median + LPF)');
    xlabel('Zaman (s)'); ylabel('Amplitude'); axis tight; grid on;
    
    saveas(gcf, fullfile(output_folder, 'ANALYSIS_IMPULSIVE.png'));
end

function seg_snr = calc_seg_snr(clean, denoised, fs)
    % Sinyali 20ms'lik segmentlere bölerek SNR hesaplar
    frame_len = round(0.02 * fs);
    num_frames = floor(length(clean) / frame_len);
    snr_vals = zeros(num_frames, 1);
    
    for i = 1:num_frames
        idx = (i-1)*frame_len + 1 : i*frame_len;
        sig_pwr = sum(clean(idx).^2);
        err_pwr = sum((clean(idx) - denoised(idx)).^2);
        snr_vals(i) = 10 * log10(sig_pwr / (err_pwr + eps));
    end
    % Aşırı düşük değerleri (sessizlik bölgeleri) temizle ve ortala
    snr_vals = snr_vals(snr_vals > -10 & snr_vals < 35);
    seg_snr = mean(snr_vals);
end