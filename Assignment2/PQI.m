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

wrist_ppg   = record(2,:);
f_sample = hdr.frequency(2);

pQI = PQI_func(wrist_ppg, f_sample);

function pulseQualityIndex = PQI_func(rawPPGsignal, samplingRate)
    
    % Call Pre-Processing
    [PPG_slimBand, PPGwideBand] = preprocessing(rawPPGsignal, samplingRate);
    test
    
%     % call (...)
%     [PpTemp, amplitudeCorrectionFactors, PpPQI, BeatTimes] = PPGsegmentationAndBeatLocalization(PPG_slimBand, PPGwideBand);
%     
%     % call (...)
%     % per each pulse -> loop?
%     templateAD = templateCreartion(PpTemp, amplitudeCorrectionFactors, PpPQI);
%     
%     % call (...)
%     pulseQualityIndex = pulseTemplateComparision(templateAD);
    
end
    
function [PPG_slimBand, PPGwideBand] = preprocessing(rawPPGsignal, samplingRate)

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

    time = ((1:size(record, 2)))/f_sample; % needed?

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
    plot(time, PPG_slimBand);
    title('wideBP filtered');
    xlabel('time [s]');
    ylabel('Amplitude [?]');
    hold off;

end

function [PpTemp, amplitudeCorrectionFactors, PpPQI, BeatTimes] = PPGsegmentationAndBeatLocalization(PPG_slimBand, PPGwideBand)

%     PpTemp = 
%     amplitudeCorrectionFactors =
%     PpPQI =
%     BeatTimes = % time of detected beats (in seconds after the start of the signal

end
