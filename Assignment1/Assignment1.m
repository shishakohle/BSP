clear;
close all;
clc;

%% Comment and uncomment according to your EDF file location

% choose on an EDF subfolder:
% edfSubfolder = ".";
edfSubfolder = "edf";

% choose on an EDF file:
edfFilename  = "s1_high_resistance_bike.edf";
% edfFilename  = "";
% edfFilename  = "";

% choose on the path notion of your operating system:
%filepath = edfSubfolder + "\" + edfFilename; % use this line for MS Windows
filepath = edfSubfolder + "/" + edfFilename; % use this line for Linux distributions / Mac OS X

%% Execute the essence function of Assignment 1

% features = extract_basic_features(filepath);


%% Function extract_basic_features as a solution to Assignment 1
% argument: TODO
% return value(s): TODO

% function features = extract_basic_features(filepath)

    % Read header and rawdata from the EDF file
    [hdr, record] = edfread(filepath);
    
    wrist_ppg   = record(2,:);
    f_sample = hdr.frequency(2);
    
    % plot original signal
    % figure;
    % time = ((1:size(record, 2)) - 0.5) / sample_freq;
    % plot(time, wrist_ppg);
    
    %% Filter signal
    
    % apply bandpass
    min_freq =  30/60; % minimum expectable frequency is 30 bpm -> converted to Hz
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    % these values are assumed to cover the range for every human individual...
    % (healthy/unhelthy/extrtaordinarily fit, young, old etc.)
    % they could be narrowed down according to a certain use case
    
    % wrist_ppg_filtered = bandpass(wrist_ppg, [min_freq max_freq], sample_freq);
    
    fpass = [min_freq max_freq]/(f_sample*0.5); 
    [b, a] = butter(2, fpass, 'bandpass'); % 2nd order BP seems to serve the purpose best
    wrist_ppg_filtered = filter(b, a, wrist_ppg);
    time = ((1:size(record, 2)))/f_sample;
    
    % optionally plot unfiltered and filtered signal
%     figure;
%     hold on;
%     subplot(2, 1, 1);
%     plot(time, wrist_ppg);
%     title('Unfiltered');
%     xlabel('time [s]');
%     ylabel('Amplitude [?]');
%     subplot(2, 1, 2);
%     plot(time, wrist_ppg_filtered);
%     title('BP filtered');
%     xlabel('time [s]');
%     ylabel('Amplitude [?]');
%     hold off;
    
    %% Segment signal into 60s windows
    
    % calculate some characteristics of the segments, that will be used
    % throughout the following on a regular basis
    width_wrist_ppg_filtered = length(wrist_ppg_filtered);
    duration_window = 60; % duration of a window in seconds
    datapoints_per_window = f_sample * duration_window;
    window_count = ceil( width_wrist_ppg_filtered / datapoints_per_window );
    
    % segment the PPG signal into 60s windows
    wrist_ppg_filtered_windows = transpose( buffer(wrist_ppg_filtered, datapoints_per_window, 0) );
    
    % replace trailing zeros in the last window with NaN, using modulo
    % (otherwise mean and var would be wrong)
    wrist_ppg_filtered_windows(end, mod(length(wrist_ppg_filtered),datapoints_per_window)+1:end) = NaN;
    
    %% Extract features
    
    % preallocating matrix for the features with NaN values - one row for each time window and 8 columns for 8 extracted features each
    features = NaN(window_count, 8);
    
    %% statistical features (time domain) (Assignment 1.1)
    features(:, 1) = mean(wrist_ppg_filtered_windows, 2,    'omitnan'); % mean of all rows, ignoring NaN values
    features(:, 2) = var (wrist_ppg_filtered_windows, 0, 2, 'omitnan'); % variance of all rows, with default weighting, ignoring NaN values
    
    %% frequency domain features based on fourier transform (Assignment 1.2)
    
    % provisory choose a single signal
    signal = wrist_ppg_filtered_windows(8,:);
    
    % calculate Power Spectral Density (PSD), following the hints in:
    % https://de.mathworks.com/help/matlab/math/fft-for-spectral-analysis.html
    % and
    % https://stackoverrun.com/de/q/4139191
    
    % calculte DFT (i.e. the complex Fourier Coefficients)
    % by utilizing the FFT algorithm
    DFT = fft(signal);
    
    % "Compute the power spectral density, a measurement of the energy at
    %  various frequencies, using the complex conjugate (CONJ)."
    PSD = DFT .* conj(DFT) / length(signal); % TODO remove NaNs from signal first!!
    
    % "Form a frequency axis for the first [...] points and use it to plot
    %  the result. (The remainder of the points are symmetric.)"
    f_axis = f_sample / length(signal) * (0:length(signal)/2);
    % f_axis and PSD: whats their correct length? they must be the same
    % length for sure. how can we properly assign a frequency to every
    % point in the PDS vector?
    
    % plot the PSD
    figure;
    plot(f_axis, PSD(1:length(signal)/2+1));
    % f_axis and PSD: whats their correct length? they must be the same
    % length for sure.
    title('Power spectral density');
    xlabel('Frequency (Hz)');
    
    %% find the frequency with the maximum power (Assignment 1.2a)
    [powerMax, f_max_index] = max(PSD);
    f_max_a = f_axis(f_max_index+1);
    % f_max_a seems to contain f_max_b, but shifted +1 in index??
    f_max_b = f_max_index * f_sample / length(signal);
    disp("f_max_a = " + f_max_a + " Hz. f_max_b = " + f_max_b + " Hz.");
    
    % add the frequency with the max. pwer to the features
    % TODO
    
    %% find the median frequency (Assignment 1.2b)
    
    % Ingo's approach
    % TODO
    
    % add the median frequency to features
    % TODO
    
    %% find the mean frequency (Assignment 1.2c)
    
    % "freq = meanfreq(pxx,f) returns the mean frequency of a power
    %  spectral density (PSD) estimate, pxx. The frequencies, f,
    %  correspond to the estimates in pxx."
    % TODO
    
    % Ingo's approach
    % TODO
    
    % add the mean frequency to features
    % TODO
    
    %% time domain features based on inter-beat intervals (Assignment 1.3)
    min_RRinterval = 1/max_freq; %minimum possible time period between heartbeats [s]

    peak_intervals = NaN(window_count, datapoints_per_window); % preallocating matrix for 60 second windows with NaN values

    for i = 1 : size(wrist_ppg_filtered_windows,1)
        signal = wrist_ppg_filtered_windows(i, :);
        time   = time( 1:length(signal) );
        [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(signal, time, 'MinPeakDistance', min_RRinterval);
        length_peak_widths = size(peak_widths, 2);
        peak_intervals(i, 1:length_peak_widths) = peak_widths; %peak intervals [s]
        
        % optionally plot signal windows with peaks
%         figure;
%         findpeaks(signal, time, 'MinPeakDistance', min_RRinterval); %finds peaks in the signal which have to be seperated by at least the max expectable frequency
%         title('PPG signal - peaks')
%         xlabel('Time [s]')
%         ylabel('Amplitude [?]')
    end
    
    features(:, 6) = mean(peak_intervals, 2, 'omitnan');
    features(:, 7) = var(peak_intervals, 0, 2, 'omitnan');
    features(:, 8) = skewness(peak_intervals, 0, 2);
    
%end
