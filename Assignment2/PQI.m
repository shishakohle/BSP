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
filepath = edfSubfolder + "\" + edfFilename;
% uncomment the next line for Linux distributions / Mac OS X
% filepath = edfSubfolder + "/" + edfFilename;

clear edfSubfolder;
clear edfFilename;

%% aquire signal
[hdr, record] = edfread(filepath);

% besser implementieren - wie in feedback damals
rawPPGsignal   = record(2,:);
samplingRate = hdr.frequency(2);


pqi = PQI_func(rawPPGsignal, samplingRate, record);

function pulseQualityIndex = PQI_func(rawPPGsignal, samplingRate, record)
    
%     Call Pre-Processing
    [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, samplingRate);
    
    % call (...)
    [PpTemp, amplitudeCorrectionFactors, PpPQI, BeatTimes] = PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, samplingRate, record);
    
    % call (...)
    % per each pulse -> loop?
    templateAD = templateCreartion(PpTemp, amplitudeCorrectionFactors, PpPQI);
    
    % call (...)
    pulseQualityIndex = pulseTemplateComparision(templateAD);
    
end
    
function [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, samplingRate)

    %% Filter signal % acc. to Papini et al 2018
    
    % apply bandpass
    min_freq =  0.4;
    max_freq_slimBand = 2.25;
    max_freq_wideBand = 10;
    
    fpass_slimBand = [min_freq max_freq_slimBand]/(samplingRate*0.5);
    fpass_wideBand = [min_freq max_freq_wideBand]/(samplingRate*0.5);
    [b_slimBand, a_slimBand] = butter(3, fpass_slimBand, 'bandpass');
    [b_wideBand, a_wideBand] = butter(3, fpass_wideBand, 'bandpass');
    PPG_slimBand = filter(b_slimBand, a_slimBand, rawPPGsignal);
    PPG_wideBand = filter(b_wideBand, a_wideBand, rawPPGsignal);
    
    %% optionally plot each signal
    
    % choose on an EDF subfolder:
    % edfSubfolder = ".";
    edfSubfolder = "edf";

    % choose on an EDF file:
    edfFilename  = "s1_high_resistance_bike.edf";
    % edfFilename  = ""; % TODO add further files
    % edfFilename  = "";

    % choose on the path notion of your operating system:
    % uncomment the next line for MS Windows
    filepath = edfSubfolder + "\" + edfFilename;
    % uncomment the next line for Linux distributions / Mac OS X
    % filepath = edfSubfolder + "/" + edfFilename;

    clear edfSubfolder;
    clear edfFilename;
    
    [hdr, record] = edfread(filepath);

    wrist_ppg   = record(2,:);
    f_sample = hdr.frequency(2);

    time = ((1:size(record, 2)))/f_sample;

    figure;
    hold on;
    subplot(3, 1, 1);
    plot(time, wrist_ppg);
    title('Unfiltered');
    xlabel('time [s]');
    ylabel('Amplitude [?]');
    subplot(3, 1, 2);
    plot(time, PPG_slimBand);
    title('slimBP filtered');
    xlabel('time [s]');
    ylabel('Amplitude [?]');
    subplot(3, 1, 3);
    plot(time, PPG_wideBand);
    title('wideBP filtered');
    xlabel('time [s]');
    ylabel('Amplitude [?]');
    hold off;
    
    figure;
    hold on;
    plot(time, PPG_slimBand);
    plot(time, PPG_wideBand);

end

function [PpTemp, amplitudeCorrectionFactors, PpPQI, BeatTimes] = PPGsegmentationAndBeatLocalization(PPG_slimBand, PPGwideBand, samplingRate, record)

%     PpTemp = 
%     amplitudeCorrectionFactors =
%     PpPQI =

    % minimum possible time period between heartbeats [s]
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    min_RRinterval = 1/max_freq;
    
    time = ((1:size(record, 2)))/samplingRate;
    time   = time( 1:length(PPG_slimBand) );
    
    [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(PPG_slimBand, time, 'MinPeakDistance', min_RRinterval);
    length_peak_widths = size(peak_widths, 2);
        
%     BeatTimes = % time of detected beats (in seconds after the start of the signal

end
