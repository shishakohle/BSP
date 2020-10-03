clc;
clear;

%% Comment and uncomment according to your EDF file location

% choose on an EDF subfolder:
% edfSubfolder = ".";
edfSubfolder = "edf";

% choose on an EDF file:
edfFilename  = "s1_high_resistance_bike.edf";
% edfFilename  = "";
% edfFilename  = "";

% choose on the path notion of your operating system:
% filepath = edfSubfolder + "\" + edfFilename; % use this line for MS Windows
filepath = edfSubfolder + "/" + edfFilename; % use this line for Linux distributions / Mac OS X

%% Execute the essence function of Assignment 1

features = extract_basic_features(filepath);


%% Function extract_basic_features as a solution to Assignment 1
% argument: TODO
% return value(s): TODO

function features = extract_basic_features(filepath)

    % Read header and rawdata from the EDF file
    [hdr, record] = edfread(filepath);
    
    wrist_ppg = record(2,:);
    sample_freq = hdr.frequency(2);
    
    % plot original signal
    figure;
    time = ((1:size(record, 2)) - 0.5) / sample_freq;
    plot(time, wrist_ppg);
    
    %% Filter and segment signal
    % apply bandpass
    min_freq = 30/60; % minimum expectable frequency is 30 bpm -> converted to Hz
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    fpass = [min_freq max_freq]; 
    
    wrist_ppg_filtered = bandpass(wrist_ppg, fpass, sample_freq);
    
    % plot filtered signal
    figure;
    time = ((1:size(record, 2))) / sample_freq;
    plot(time, wrist_ppg_filtered);
    
    % segment into 60s windows
    datapoints_per60s = sample_freq*60;
    width_wrist_ppg_filtered = size(wrist_ppg_filtered, 2);
    rownumber = size((1 : datapoints_per60s : width_wrist_ppg_filtered), 2);
    length_lastrow = mod(width_wrist_ppg_filtered, datapoints_per60s);

    wrist_ppg_filtered_windowed = NaN(rownumber, datapoints_per60s); % preallocating matrix for 60 second windows with NaN values
    
    count = 1;
    for i = 1 : datapoints_per60s : width_wrist_ppg_filtered
        if (i-1+datapoints_per60s) > width_wrist_ppg_filtered
            wrist_ppg_filtered_windowed(count, 1:length_lastrow) = wrist_ppg_filtered(1, i:end);
            break;
        end
        wrist_ppg_filtered_windowed(count, :) = wrist_ppg_filtered(1, i:(i-1+datapoints_per60s));
        count = count+1;
    end
    
    %% Extract features
    features = NaN(rownumber, 8); % preallocating matrix for the features with NaN values - one row per time window and 8 columns for the extracted features
    
    % statistical features time domain
    features(:, 1) = mean(wrist_ppg_filtered_windowed, 2, 'omitnan'); % mean of all rows, ignoring NaN values
    features(:, 2) = var(wrist_ppg_filtered_windowed, 0, 2, 'omitnan'); % variance of all rows, with default weighting, ignoring NaN values
    
    % frequency domain features based on fourier transform
        % use fft() function ?
    
    % time domain features based on inter-beat intervals
    
end