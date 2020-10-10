%% Biosignal Processing - Assingment 1
% Group H: Alexander Neumayer, Lukas Riedler, Laura Kummer, Ingo Weigel
% 2020-10-10
% Find this code in the project repository on GitHub:
% https://github.com/shishakohle/BSP

clear;
close all;
clc;

%% Comment and uncomment according to your EDF file location and your OS

% choose on an EDF subfolder:
% edfSubfolder = ".";
edfSubfolder = "edf";

% choose on an EDF file:
edfFilename  = "s1_high_resistance_bike.edf";
% edfFilename  = ""; % TODO add further files
% edfFilename  = "";

% choose on the path notion of your operating system:
% uncomment the next line for MS Windows
% filepath = edfSubfolder + "\" + edfFilename;
% uncomment the next line for Linux distributions / Mac OS X
filepath = edfSubfolder + "/" + edfFilename;

clear edfSubfolder;
clear edfFilename;

%% Execute the essence function of Assignment 1

features = extract_basic_features(filepath);


%% Function extract_basic_features as a solution to Assignment 1
% argument:
%   Path to an EDF file containing a Wrist PPG signal in index 2.
% return value(s):
%   A matrix containing one line per 60 seconds in the PPG signal.
%   Each column represents a certain feature of the corresponding
%   60 seconds snippet:
%       column 1: mean
%       column 2: variance
%       column 3: frequency with the maximum power
%       column 4: median frequency
%       column 5: mean frequency
%       column 6: mean inter-beat interval
%       column 7: variance of inter-beat intervals
%       column 8: skewness of the distribution of inter-beat intervals

function features = extract_basic_features(filepath)

    %% Read header and rawdata from the EDF file
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
    
    % The matrix features is to be used as the final return value in the
    % end.
    
    % preallocating matrix for the features with NaN values - one row for
    % each time window and 8 columns for 8 extracted features each
    features = NaN(window_count, 8);
    
    %% statistical features (time domain) (Assignment 1.1)
    % mean of all rows, ignoring NaN values
    features(:, 1) = mean(wrist_ppg_filtered_windows, 2,    'omitnan');
    % variance of all rows, with default weighting, ignoring NaN values
    features(:, 2) = var (wrist_ppg_filtered_windows, 0, 2, 'omitnan');
    
    %% frequency domain features based on fourier transform (Assignment 1.2)
    
    % prepare signal in such a way that it can be processed by fft() and others
    signal = transpose (wrist_ppg_filtered_windows);
    
    % calculate Power Spectral Density (PSD), following the hints in:
    % https://de.mathworks.com/help/matlab/math/fft-for-spectral-analysis.html
    % and
    % https://stackoverrun.com/de/q/4139191
    
    % calculte DFT (i.e. the complex Fourier Coefficients)
    % by utilizing the FFT algorithm
    DFT = fft(signal);
    % TODO remove NaNs from last window first!! not done yet.
    
    % "Compute the power spectral density, a measurement of the energy at
    %  various frequencies, using the complex conjugate (CONJ)."
    PSD = DFT .* conj(DFT) / length(signal);
    
    % "Form a frequency axis (...)"
    f_axis = f_sample / length(signal) * (0:length(signal)/2);
    % f_axis and PSD: whats their correct lengths?
    
    %% find the frequency with the maximum power (Assignment 1.2a)
    [powerMax, f_max_index] = max(PSD);
    f_max_a = f_axis(f_max_index+1);
    % f_max_a seems to contain f_max_b, but shifted +1 in index??
    f_max_b = f_max_index * f_sample / length(signal);
    consoleOutputs = "f_max_a = " + f_max_a + " Hz. f_max_b = " + f_max_b + " Hz.";
    disp(consoleOutputs(:));
    
    % add the frequency with the maximum power to the features
    features(:,3) = f_max_a;
    
    %% find the median frequency (Assignment 1.2b)
    f_med = medfreq(wrist_ppg_filtered_windows(:,1:9),f_sample);
    f_med(10) = medfreq(wrist_ppg_filtered_windows(~isnan(wrist_ppg_filtered_windows(:,10))),f_sample);
    % TODO this is hard-coded for the signal having 10 windows, where the
    % last window (no. 10) has NaNs to be removed from it.
    
    % add the median frequency to features
    features(:,4) = f_med;
    
    %% find the mean frequency (Assignment 1.2c)
    f_mean = meanfreq(wrist_ppg_filtered_windows(:,1:9),f_sample);
    f_mean(10) = meanfreq(wrist_ppg_filtered_windows(~isnan(wrist_ppg_filtered_windows(:,10))),f_sample);
    % TODO this is hard-coded for the signal having 10 windows, where the
    % last window (no. 10) has NaNs to be removed from it.
    
    % add the mean frequency to features
    features(:,5) = f_mean;
    
    %% time domain features based on inter-beat intervals (Assignment 1.3)
    
    % minimum possible time period between heartbeats [s]
    min_RRinterval = 1/max_freq;

    % preallocating matrix for 60 second windows with NaN values
    peak_intervals = NaN(window_count, datapoints_per_window);

    for i = 1 : size(wrist_ppg_filtered_windows,1)
        signal = wrist_ppg_filtered_windows(i, :);
        time   = time( 1:length(signal) );
        [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(signal, time, 'MinPeakDistance', min_RRinterval);
        length_peak_widths = size(peak_widths, 2);
        % peak intervals [s]
        peak_intervals(i, 1:length_peak_widths) = peak_widths;
        
        % optionally plot signal windows with peaks
%         figure;
%         findpeaks(signal, time, 'MinPeakDistance', min_RRinterval); %finds peaks in the signal which have to be seperated by at least the max expectable frequency
%         title('PPG signal - peaks')
%         xlabel('Time [s]')
%         ylabel('Amplitude [?]')
    end
    
    % add features of the inter-beat intervalls to features
    features(:, 6) = mean(peak_intervals, 2, 'omitnan');
    features(:, 7) = var(peak_intervals, 0, 2, 'omitnan');
    features(:, 8) = skewness(peak_intervals, 0, 2);
    
end
