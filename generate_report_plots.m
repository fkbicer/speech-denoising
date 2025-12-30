function generate_report_plots(s_clean, s_noisy, s_denoised, fs, noise_name)
    % s_clean: Orijinal temiz sinyal
    % s_noisy: Gürültülü sinyal
    % s_denoised: Filtrelediğin sinyal
    % fs: Sample rate (16000)
    % noise_name: Gürültü türü (Örn: 'Babble', 'White')

    % SNR Hesaplamaları (Ödev Task 5)
    snr_pre = 10 * log10(sum(s_clean.^2) / sum((s_clean - s_noisy).^2));
    snr_post = 10 * log10(sum(s_clean.^2) / sum((s_clean - s_denoised).^2));
    gain = snr_post - snr_pre;

    % Figür Oluşturma
    fig = figure('Name', ['Report Analysis - ' noise_name], 'Color', 'w', 'Position', [50 50 1200 900]);
    t = (0:length(s_clean)-1)/fs;

    % --- TASK 2 & 6: TIME DOMAIN WAVEFORMS ---
    subplot(3, 2, 1);
    plot(t, s_noisy, 'Color', [0.7 0.7 0.7]); hold on;
    plot(t, s_clean, 'b', 'LineWidth', 0.8);
    title([noise_name ' - Time Domain (Noisy vs Clean)']);
    ylabel('Amplitude'); xlabel('Time (s)'); legend('Noisy', 'Clean');
    grid on; xlim([1 1.5]); % Zoomed for visibility

    subplot(3, 2, 2);
    plot(t, s_denoised, 'r', 'LineWidth', 0.8); hold on;
    plot(t, s_clean, 'b', 'LineWidth', 0.5);
    title(['Denoised (SNR Post: ' num2str(snr_post, '%.2f') ' dB)']);
    ylabel('Amplitude'); xlabel('Time (s)'); legend('Denoised', 'Clean');
    grid on; xlim([1 1.5]);

    % --- TASK 2: FREQUENCY DOMAIN (Magnitude Spectra) ---
    [f_vec, P_noisy] = get_spectrum(s_noisy, fs);
    [~, P_denoised] = get_spectrum(s_denoised, fs);
    [~, P_clean] = get_spectrum(s_clean, fs);

    subplot(3, 2, 3);
    plot(f_vec/1000, P_noisy, 'Color', [0.7 0.7 0.7]); hold on;
    plot(f_vec/1000, P_clean, 'b');
    title('Frequency Spectrum (Noisy)');
    ylabel('|X(e^{j\omega})| (dB)'); xlabel('Frequency (kHz)');
    grid on; xlim([0 fs/2000]);

    subplot(3, 2, 4);
    plot(f_vec/1000, P_denoised, 'r'); hold on;
    plot(f_vec/1000, P_clean, 'b');
    title(['Frequency Spectrum (Denoised) - Gain: ' num2str(gain, '%.2f') ' dB']);
    ylabel('|X(e^{j\omega})| (dB)'); xlabel('Frequency (kHz)');
    grid on; xlim([0 fs/2000]);

    % --- TASK 6: SPECTROGRAMS ---
    subplot(3, 2, 5);
    draw_spec(s_noisy, fs, [noise_name ' Noisy Spectrogram']);
    
    subplot(3, 2, 6);
    draw_spec(s_denoised, fs, [noise_name ' Denoised Spectrogram']);

    % Alt bilgi ekleme
    sgtitle(['DSP Speech Denoising Project: ' noise_name ' Analysis'], 'FontSize', 16, 'FontWeight', 'bold');
end

% Yardımcı Fonksiyon: Spektrum Hesaplama
function [f, P1] = get_spectrum(x, fs)
    L = length(x);
    Y = fft(x);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(floor(L/2)))/L;
    P1 = 20*log10(P1 + eps); % dB scale
end

% Yardımcı Fonksiyon: Spektrogram Çizimi
function draw_spec(x, fs, t_title)
    [s, f, t] = spectrogram(x, hamming(512), 256, 512, fs, 'yaxis');
    imagesc(t, f/1000, 20*log10(abs(s)+eps));
    axis xy; colormap jet; caxis([-90 0]);
    title(t_title); ylabel('Freq (kHz)'); xlabel('Time (s)');
    colorbar;
end