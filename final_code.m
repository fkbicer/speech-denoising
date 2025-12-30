clc; clear; close all;

% 1. SETUP
fs_target = 16000;
file_clean = 'clean_speech.wav';

% Filename and Noise Type mapping
file_list = {
    'clean_white_5dB.wav',     'WHITE';
    'clean_street_5dB.wav',    'STREET';
    'clean_pink_5dB.wav',      'PINK';
    'clean_babble_5dB.wav',    'BABBLE';
    'clean_impulsive_5dB.wav', 'IMPULSIVE'
};

% Check if reference file exists
if ~isfile(file_clean), error('Clean speech file missing!'); end

[x_clean, fs] = audioread(file_clean);
% Force 16k sample rate
if fs ~= fs_target, x_clean = resample(x_clean, fs_target, fs); fs=fs_target; end

fprintf('=== DSP PROJECT STARTED ===\n');

% 2. PROCESSING LOOP
for i = 1:size(file_list, 1)
    target_file = file_list{i, 1};
    noise_type  = file_list{i, 2};
    
    if ~isfile(target_file)
        fprintf('Skipping: %s (Not found)\n', target_file);
        continue;
    end
    
    fprintf('\n>>> Processing: %s (%s)...\n', target_file, noise_type);
    [x_noisy, fs_in] = audioread(target_file);
    if fs_in ~= fs, x_noisy = resample(x_noisy, fs, fs_in); end
    
    % Trim signals to match lengths for SNR calc
    L = min(length(x_clean), length(x_noisy));
    x_c = x_clean(1:L);
    x_n = x_noisy(1:L);
    
    % 3. PRE-PROCESSING
    % Cut heavy rumble for Street noise immediately
    if strcmp(noise_type, 'STREET')
        [b_pre, a_pre] = butter(4, 180/(fs/2), 'high'); 
        x_proc = filtfilt(b_pre, a_pre, x_n);
    else
        x_proc = x_n;
    end
    
    % 4. FILTERING STRATEGIES
    
    if strcmp(noise_type, 'IMPULSIVE')
        % Strategy: Median Filter kills the spikes
        y_temp = medfilt1(x_proc, 3); % Window=3 is the sweet spot
        
        % Treble boost to bring back some clarity
        [b_high, a_high] = butter(2, 2500/(fs/2), 'high');
        highs = filtfilt(b_high, a_high, y_temp) * 0.30; 
        y_denoised = y_temp + highs;
        
        method_name = 'Median(3) + TrebleBoost';
        
    elseif strcmp(noise_type, 'PINK')
        % Strategy: Balance the spectrum then smooth
        
        % 1. High-pass to fix the 1/f bass heaviness
        [b_hp, a_hp] = butter(4, 150/(fs/2), 'high');
        x_pink = filtfilt(b_hp, a_hp, x_proc);
        
        % 2. Wiener Filter in freq domain
        win = 512; hop = 256; nfft = 512;
        [S, ~, ~] = stft(x_pink, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        mag = abs(S); phs = angle(S);
        
        mag_denoised = wiener2(mag, [7 7]);
        
        % Simple masking
        mask = mag_denoised > (0.05 * max(mag_denoised(:)));
        mag_denoised = mag_denoised .* mask;
        
        % Reconstruct
        S_new = mag_denoised .* exp(1j * phs);
        y_temp = istft(S_new, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        y_temp = real(y_temp);
        
        % 3. Air Boost to keep it natural
        [b_air, a_air] = butter(2, 4000/(fs/2), 'high');
        y_denoised = y_temp + (filtfilt(b_air, a_air, y_temp) * 0.15);
        
        method_name = 'HPF(150Hz) + Wiener + Air';

    elseif strcmp(noise_type, 'WHITE')
        % Strategy: Cut the hiss
        
        % 1. LPF at 5.5kHz 
        [b_lp, a_lp] = butter(6, 5500/(fs/2), 'low'); 
        x_white = filtfilt(b_lp, a_lp, x_proc);
        
        % 2. Standard Wiener
        win = 512; hop = 256; nfft = 512;
        [S, ~, ~] = stft(x_white, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        mag = abs(S); phs = angle(S);
        
        mag_denoised = wiener2(mag, [5 5]);
        mask = mag_denoised > (0.06 * max(mag_denoised(:)));
        mag_denoised = mag_denoised .* mask;
        
        S_new = mag_denoised .* exp(1j * phs);
        y_temp = istft(S_new, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        y_denoised = real(y_temp);
        
        % 3. Air Boost
        [b_air, a_air] = butter(2, 4000/(fs/2), 'high');
        y_denoised = y_denoised + (filtfilt(b_air, a_air, y_denoised) * 0.15);
        
        method_name = 'LPF(5.5k) + Wiener + Air';
        
    else
        % Strategy: Standard Wiener for Street/Babble
        
        win = 512; hop = 256; nfft = 512;
        [S, ~, ~] = stft(x_proc, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        mag = abs(S); phs = angle(S);
        
        mag_denoised = wiener2(mag, [7 7]);
        
        % Harder gate for Babble since speech overlaps
        if strcmp(noise_type, 'BABBLE'), thresh = 0.08; else, thresh = 0.04; end
        
        mask = mag_denoised > (thresh * max(mag_denoised(:)));
        mag_denoised = mag_denoised .* mask;
        
        S_new = mag_denoised .* exp(1j * phs);
        y_temp = istft(S_new, fs, 'Window', hamming(win), 'OverlapLength', hop, 'FFTLength', nfft);
        y_denoised = real(y_temp);
        
        method_name = 'Wiener + Gate';
    end
    
    % 5. POST-PROCESSING
    
    % Fix array size mismatch
    if length(y_denoised) > L, y_denoised = y_denoised(1:L);
    elseif length(y_denoised) < L, y_denoised = [y_denoised; zeros(L-length(y_denoised),1)]; end
    
    % General cleanup for low freq rumble (except Impulsive/Pink)
    if ~strcmp(noise_type, 'IMPULSIVE') && ~strcmp(noise_type, 'PINK')
        [b_hp, a_hp] = butter(4, 120/(fs/2), 'high');
        y_denoised = filtfilt(b_hp, a_hp, y_denoised);
    end
    
    % Focus on speech band for Babble
    if strcmp(noise_type, 'BABBLE')
        [b_bp, a_bp] = butter(6, [300 3400]/(fs/2), 'bandpass');
        y_denoised = filtfilt(b_bp, a_bp, y_denoised);
    end
    
    % Normalize volume so it doesn't clip
    m = max(abs(y_denoised));
    if m > 0, y_denoised = y_denoised / m * 0.95; end
    
    % 6. METRICS & SAVING
    snr_pre = 10 * log10(sum(x_c.^2) / sum((x_c - x_n).^2));
    snr_post = 10 * log10(sum(x_c.^2) / sum((x_c - y_denoised).^2));
    gain = snr_post - snr_pre;
    
    fprintf('   Method: %s\n', method_name);
    fprintf('   SNR: %.2f dB -> %.2f dB (Gain: +%.2f dB)\n', snr_pre, snr_post, gain);
    
    out_name = sprintf('DENOISED_%s.wav', target_file(1:end-4));
    audiowrite(out_name, y_denoised, fs);
    fprintf('   [AUDIO] Saved: %s\n', out_name);
    
    plot_comparative(x_n, y_denoised, fs, noise_type, snr_pre, snr_post);
    fig_name = sprintf('ANALYSIS_%s.png', noise_type);
    saveas(gcf, fig_name);
    fprintf('   [PLOT] Saved: %s\n', fig_name);
    
    drawnow;
end

% 7. REFERENCE PLOT
fprintf('\n>>> Plotting Clean Reference...\n');
plot_clean_reference(x_clean, fs);
saveas(gcf, 'REFERENCE_CLEAN.png');
fprintf('   [REF] Saved.\n');
fprintf('\n=== DONE ===\n');

%  FUNCTIONS

function plot_comparative(noisy, denoised, fs, type_name, snr_in, snr_out)
    figure('Name', [type_name ' Analysis'], 'Color', 'w', 'Position', [50 50 1000 800]);
    t = (0:length(noisy)-1)/fs;
    
    % Waveform
    subplot(3,2,1); plot(t, noisy, 'Color', [0.8 0.2 0.2]); axis tight; ylim([-1 1]); grid on;
    title([type_name ' - Noisy (SNR=' num2str(snr_in,'%.1f') 'dB)'], 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Amplitude');
    
    subplot(3,2,2); plot(t, denoised, 'Color', [0.2 0.6 0.2]); axis tight; ylim([-1 1]); grid on;
    title([type_name ' - Denoised (SNR=' num2str(snr_out,'%.1f') 'dB)'], 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Amplitude');
    
    % Spectrum
    [f_kHz, P_n] = calc_raw_fft_db(noisy, fs);
    [~, P_c] = calc_raw_fft_db(denoised, fs);
    
    subplot(3,2,3); plot(f_kHz, P_n, 'Color', [0.8 0.2 0.2]); xlim([0 8]); grid on;
    title('Noisy Spectrum', 'FontSize', 10, 'FontWeight', 'bold'); ylabel('Mag (dB)');
    
    subplot(3,2,4); plot(f_kHz, P_c, 'Color', [0.2 0.6 0.2]); xlim([0 8]); grid on;
    title('Denoised Spectrum', 'FontSize', 10, 'FontWeight', 'bold'); ylabel('Mag (dB)');
    
    % Spectrograms
    subplot(3,2,5); draw_spectrogram(noisy, fs, 'Noisy Spectrogram');
    subplot(3,2,6); draw_spectrogram(denoised, fs, 'Denoised Spectrogram');
end

function plot_clean_reference(clean, fs)
    figure('Name', 'Clean Reference', 'Color', 'w', 'Position', [150 150 800 900]);
    t = (0:length(clean)-1)/fs;
    
    subplot(3,1,1); plot(t, clean, 'b'); axis tight; ylim([-1 1]); grid on;
    title('REFERENCE: Clean Speech Waveform', 'FontSize', 12, 'FontWeight', 'bold');
    
    [f_kHz, P_dB] = calc_raw_fft_db(clean, fs);
    subplot(3,1,2); plot(f_kHz, P_dB, 'b'); xlim([0 8]); grid on;
    title('REFERENCE: Clean Speech Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
    
    subplot(3,1,3); draw_spectrogram(clean, fs, 'REFERENCE: Clean Speech Spectrogram');
end

function [f_kHz, P_dB] = calc_raw_fft_db(x, fs)
    L = length(x); Y = fft(x);
    P2 = abs(Y); 
    P1 = P2(1:floor(L/2)+1); P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(floor(L/2)))/L;
    f_kHz = f/1000;
    P_dB = 20*log10(P1 + eps);
end

function draw_spectrogram(x, fs, t_title)
    % High-Res settings (looks cleaner)
    win_size = 1024; 
    noverlap = 512; 
    nfft = 1024;
    
    [s, f, t] = spectrogram(x, hamming(win_size), noverlap, nfft, fs);
    S_dB = 20*log10(abs(s));
    
    % Normalize peak to 0dB
    max_val = max(max(S_dB));
    S_norm = S_dB - max_val;
    
    imagesc(t, f/1000, S_norm); 
    axis xy; 
    colormap parula; 
    
    c = colorbar;
    c.Label.String = 'Norm. Power (dB)';
    
    % Contrast fix: -80dB floor gives a nice deep blue background
    caxis([-80 0]); 
    
    title(t_title, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time (s)'); ylabel('Freq (kHz)');
end