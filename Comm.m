%% Parameters

Fs_standard = 44100; % Standard sampling rate
Fc1 = 100e3; % Carrier frequency for signal1 (100 kHz)
Fc2 = 150e3; % Carrier frequency for signal2 (150 kHz)
upsampling_factor = 10; % Factor to upsample the signal
IF_freq = 25e3; % Intermediate Frequency (25 kHz)
F_shift1 = 0.2e3; % Small frequency shift (0.2 kHz)
F_shift2 = 1.2e3; % Large frequency shift (1.2 kHz)

%% Transmitter

% Read and Analyze the Input Signals
[signal1, Fs] = audioread('Short_FM9090.wav');
mono_signal1 = sum(signal1, 2);

[signal2, Fs] = audioread('Short_BBCArabic2.wav');
mono_signal2 = sum(signal2, 2);

% zero-padding shorter signals
max_len = max([length(mono_signal1), length(mono_signal2)]);
mono_signal1 = [mono_signal1; zeros(max_len - length(mono_signal1), 1)];
mono_signal2 = [mono_signal2; zeros(max_len - length(mono_signal2), 1)];

% Analyze the spectrum using FFT
plotSpectrum(mono_signal1, Fs);
title('Spectrum of Signal 1');

plotSpectrum(mono_signal2, Fs);
title('Spectrum of Signal 2');

% Upsample the signals to avoid aliasing
mono_signal1 = interp(mono_signal1, upsampling_factor, 1);
mono_signal2 = interp(mono_signal2, upsampling_factor, 1);
Fs = Fs * upsampling_factor; % Update sampling frequency

% Analyze the spectrum using FFT
plotSpectrum(mono_signal1, Fs);
title('Spectrum of Signal 1');

plotSpectrum(mono_signal2, Fs);
title('Spectrum of Signal 2');

% Generate the carrier wave for modulation
t = (0:length(mono_signal1)-1)' / Fs; % Time vector
carrier1 = cos(2 * pi * Fc1 * t); % Carrier for signal1
carrier2 = cos(2 * pi * Fc2 * t); % Carrier for signal2

% DSB-SC Modulation
modulated_signal1 = mono_signal1 .* carrier1;
modulated_signal2 = mono_signal2 .* carrier2;

% FDM Signal
fdm_signal = modulated_signal1 + modulated_signal2;

% Plot the spectra of modulated and FDM signals
plotSpectrum(modulated_signal1, Fs);
title('Spectrum of Modulated Signal 1');

plotSpectrum(modulated_signal2, Fs);
title('Spectrum of Modulated Signal 2');

plotSpectrum(fdm_signal, Fs);
title('Spectrum of FDM Signal');

%% Receiver

% Step 1: RF Stage - Bandpass Filter around Fc1 frequency
BPF_RF = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', Fc1 - 8.5e3, 'HalfPowerFrequency2', Fc1 + 8.5e3, ...
    'SampleRate', Fs);

rf_signal1 = filter(BPF_RF, fdm_signal); % Apply bandpass filter to the FDM signal

%rf_signal1 = fdm_signal; % When removing RF stage

% Plot RF stage output spectrum
plotSpectrum(rf_signal1, Fs);
title('Spectrum of RF Stage Output');

% Step 2: Oscillator - Generate Local Oscillator at (Fc + IF)
t = (0:length(fdm_signal)-1) / Fs; % Time vector based on sampling frequency
oscillator = cos(2 * pi * (Fc1 + IF_freq) * t); % Local oscillator for mixing

rf_signal1 = rf_signal1(:); % Force column vector
oscillator = oscillator(:); % Force column vector

% Multiply RF signal with oscillator
mixed_signal = rf_signal1 .* oscillator;

% Plot mixer output spectrum
plotSpectrum(mixed_signal, Fs);
title('Spectrum of Oscillator Output');

% Step 3: IF Stage - Bandpass Filter around IF frequency
BPF_IF = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', IF_freq - 8.5e3, 'HalfPowerFrequency2', IF_freq + 8.5e3, ...
    'SampleRate', Fs);

if_signal = filter(BPF_IF, mixed_signal); % Apply bandpass filter at IF frequency

% Plot IF stage output spectrum
plotSpectrum(if_signal, Fs);
title('Spectrum of IF Stage Output');

% Step 4: Baseband Detection - Mix with Carrier and Apply Low-Pass Filter
% Generate baseband carrier for demodulation
baseband_carrier = cos(2 * pi * IF_freq * t); % Carrier at Fc1 for demodulation

if_signal = if_signal(:); % Force column vector
baseband_carrier = baseband_carrier(:); % Force column vector

% Mix the IF signal with the baseband carrier
baseband_signal = if_signal .* baseband_carrier;

% Plot Baseband mixer output spectrum
plotSpectrum(baseband_signal, Fs);
title('Spectrum of Baseband Signal Before LPF');

% Design Low-Pass Filter for baseband signal
LPF_baseband = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 8.5e3, 'SampleRate', Fs);

% Apply LPF to extract baseband signal
baseband_output = filter(LPF_baseband, baseband_signal); 

% Resample the signal to standard audio rate
baseband_output_resampled = resample(baseband_output, Fs_standard, Fs); % Downsample signal
%Fs = Fs_standard;

% Step 7: Plot and Listen to the Demodulated Signal
plotSpectrum(baseband_output, Fs);
title('Spectrum of Baseband Output After LPF');
xlim([-2.5e4 2.5e4]); % Zoom to see the signal clearly

% Play the demodulated audio
sound(baseband_output_resampled, Fs_standard); % Play the baseband signal
pause(length(baseband_output_resampled) / Fs_standard); % Wait for Signal to finish before playing another

%% Adding noise to the RF output signal

% Add AWGN noise to the demodulated signal
snr = 10; % Signal-to-Noise Ratio (SNR) in dB
noisy_rf_signal1 = awgn(rf_signal1, snr, 'measured'); % Add noise

% Plot the noisy RF signal
plotSpectrum(noisy_rf_signal1, Fs);
title('Spectrum of Noisy RF Stage Output');

% Multiply noisy RF signal with oscillator
noisy_mixed_signal = noisy_rf_signal1 .* oscillator;

noisy_if_signal = filter(BPF_IF, noisy_mixed_signal); % Apply bandpass filter at IF frequency

% Mix the noisy IF signal with the baseband carrier
noisy_baseband_signal = noisy_if_signal .* baseband_carrier;

% Apply LPF to extract noisy baseband signal
noisy_baseband_output = filter(LPF_baseband, noisy_baseband_signal); 

% Resample the noisy signal to standard audio rate
noisy_baseband_output_resampled = resample(noisy_baseband_output, Fs_standard, Fs); % Downsample signal

% Plot the noisy signal
plotSpectrum(noisy_baseband_output_resampled, Fs_standard);
title('Spectrum of Noisy Baseband Output After LPF');
xlim([-2.5e4 2.5e4]); % Zoom to see the signal clearly

% Play the noisy signal
sound(noisy_baseband_output_resampled, Fs_standard); 


%% Functions

% Function to plot the spectrum and time-domain of the signal
function plotSpectrum(signal, Fs)
    N = length(signal); % Length of the signal
    t = (0:N-1) / Fs; % Time vector
    freq = linspace(-Fs/2, Fs/2, N); % Frequency axis for the FFT
    spectrum = fftshift(abs(fft(signal))); % Compute the shifted spectrum
    
    % Create a figure with two subplots
    figure;
    
    % Subplot for time-domain signal
    subplot(2, 1, 1);
    plot(t, signal);
    title('Time-Domain Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    %ylim([-3 3]); 
    
    % Subplot for frequency-domain signal
    subplot(2, 1, 2);
    plot(freq, spectrum);
    title('Frequency-Domain Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;

    % Enable Data Cursor Mode for frequency plot
    datacursormode on; % Turns on data cursor mode for figure
end


