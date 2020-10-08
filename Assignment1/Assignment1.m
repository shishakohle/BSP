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
    
    % statistical features (time domain) (Assignment 1.1)
    features(:, 1) = mean(wrist_ppg_filtered_windows, 2,    'omitnan'); % mean of all rows, ignoring NaN values
    features(:, 2) = var (wrist_ppg_filtered_windows, 0, 2, 'omitnan'); % variance of all rows, with default weighting, ignoring NaN values
    
    % frequency domain features based on fourier transform (Assignment 1.2)
% %     signal = ppg_windows(2,:);
% %     signal = ppg_windows;
%     time   = time( 1:length(signal) );
%     signal_length = length(signal);
%     
%     figure;
%     subplot(3,1,1);
%     plot(time, signal);
%     
%     % MATLAB approach
%     % --> https://de.mathworks.com/help/matlab/math/fft-for-spectral-analysis.html
%     subplot(3,1,2);
%     % remove bias (0 Hz)??
%     % ...
%     Y = fft(signal);
%     %Pyy = Y.*conj(Y);
%     Pyy = Y.*conj(Y)/signal_length;
%     %axis_freq = sample_freq/length(signal_matlab)*(0:)
%     plot(Pyy);
%     [maxValue, indexMax] = max(Pyy);
%     f_peak = indexMax * sample_freq / signal_length;
%     disp("MATLAB approach says f_peak is " + f_peak + " Hz.");
%     %plot(freq,aproach_matlab);
%     
%     % Frequency domain
%     n = 2^nextpow2(signal_length);
%     Y_fd = fft(signal, n);
%     f = sample_freq*(0:(n/2))/n;
%     P = abs(Y/n).^2;
%     figure;
%     plot(f,P(1:n/2+1)) 
%     title('Signal in Frequency Domain')
%     xlabel('Frequency (f)')
%     ylabel('|P(f)|^2')
 
%     % FORUM approach
%     % --> https://stackoverrun.com/de/q/4139191
%     signal_forum = signal;
%     subplot(3,1,3);
%     % remove bias (0 Hz)
%     % signal_forum = signal_forum - mean(signal);
%     signal_fourier_domain = fft(signal_forum);
%     complex_magnitudes = abs(signal_fourier_domain);
%     plot(complex_magnitudes);
%     [maxValue, indexMax] = max(complex_magnitudes);
%     f_peak = indexMax * sample_freq / length(signal_forum);
%     disp("FORUM approach says f_peak is " + f_peak + " Hz.");

    % time domain features based on inter-beat intervals (Assignment 1.3)
    min_RRinterval = 1/max_freq; %minimum possible time period between heartbeats [s]

    peak_intervals = NaN(window_count, datapoints_per_window); % preallocating matrix for 60 second windows with NaN values

    for i = 1 : size(wrist_ppg_filtered_windows,1)
        signal = wrist_ppg_filtered_windows(i, :);
        time   = time( 1:length(signal) );
        [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(signal, time, 'MinPeakDistance', min_RRinterval);
        length_peak_widths = size(peak_widths, 2);
        peak_intervals(i, 1:length_peak_widths) = peak_widths; %peak intervals [s]
        
        % optionally plot signal windows with peaks
        figure;
        findpeaks(signal, time, 'MinPeakDistance', min_RRinterval); %finds peaks in the signal which have to be seperated by at least the max expectable frequency
        title('PPG signal - peaks')
        xlabel('Time [s]')
        ylabel('Amplitude [?]')
    end
    
    features(:, 6) = mean(peak_intervals, 2, 'omitnan');
    features(:, 7) = var(peak_intervals, 0, 2, 'omitnan');
    features(:, 8) = skewness(peak_intervals, 0, 2);
    
%end
