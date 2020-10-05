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
filepath = edfSubfolder + "\" + edfFilename; % use this line for MS Windows
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
    sample_freq = hdr.frequency(2);
    
    % plot original signal
    % figure;
    % time = ((1:size(record, 2)) - 0.5) / sample_freq;
    % plot(time, wrist_ppg);
    
    %% Filter and segment signal
    
    % apply bandpass
    min_freq =  30/60; % minimum expectable frequency is 30 bpm -> converted to Hz
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    
    % wrist_ppg_filtered = bandpass(wrist_ppg, [min_freq max_freq], sample_freq);
    
    fpass = [min_freq max_freq]/(sample_freq*0.5); 
    [b, a] = butter(2, fpass, 'bandpass');
    wrist_ppg_filtered = filter(b, a, wrist_ppg);
    
    % plot unfiltered and filtered signal
    figure;
    hold on;
    time = ((1:size(record, 2)))/sample_freq;
    subplot(2, 1, 1);
    plot(time, wrist_ppg);
    title('Unfiltered');
    subplot(2, 1, 2);
    plot(time, wrist_ppg_filtered);
    title('BP filtered');
    hold off;
    
    % segment into 60s windows    
    ppg_windows = transpose( buffer(wrist_ppg_filtered, sample_freq*60, 0) );
    
    %% Extract features
    features = NaN( size(ppg_windows,1), 8 ); % preallocating matrix for the features with NaN values - one row for each time window and 8 columns for 8 extracted features each
    
    % statistical features time domain (Assignment 1.1)
    features(:, 1) = mean(ppg_windows, 2,    'omitnan'); % mean of all rows, ignoring NaN values
    features(:, 2) = var (ppg_windows, 0, 2, 'omitnan'); % variance of all rows, with default weighting, ignoring NaN values
    
    % frequency domain features based on fourier transform (Assignment 1.2)
        % use fft() function ?
    
    % time domain features based on inter-beat intervals (Assignment 1.3)
    
%end
