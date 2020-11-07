function [matrix, pulseWaveTemplate] = PQI (rawPPGsignal, samplingRate, ...
    record)
    
    %% Algorithm according to Figure 1 in the Papini paper

    % call Pre-processing
    [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, ...
        samplingRate);
    
    % call PPG segmentation, beat localization
    [PP_Temp, AmplitudeCorrectionFactors, PP_PQI, BeatTimes] = ...
        PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, ...
        samplingRate, record);
    
    % call Template creation
    Temp_Ad = templateCreation(PP_Temp, AmplitudeCorrectionFactors);
    
    % call Pulse-Template comparison
    PulseQualityIndex = pulseTemplateComparision(PP_PQI, Temp_Ad);
    
    %% assemble return values according to Assignment 2 task description:
    
    % "Matrix with timing, amplitude, and quality index for each detected
    %  beat"
    % "Optional: calculated pulse wave template"
    
    matrix = [PulseQualityIndex]; % TODO: timing, amplitude
    pulseWaveTemplate = Temp_Ad; % TODO: is this accurate?
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing main blocks %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, ...
    samplingRate)
    
    % "Pre-processing" in the Papini paper, Figure 1
    
    PPG_slimBand = bandpassFilter(rawPPGsignal, 0.4, 2.25, samplingRate);
    PPG_wideBand = bandpassFilter(rawPPGsignal, 0.4, 10.0, samplingRate);
    
%     % optionally plot each signal
%     
%     % choose on an EDF subfolder:
%     % edfSubfolder = ".";
%     edfSubfolder = "edf";
% 
%     % choose on an EDF file:
%     edfFilename  = "s1_high_resistance_bike.edf";
%     % edfFilename  = ""; % TODO add further files
%     % edfFilename  = "";
% 
%     % choose on the path notion of your operating system:
%     % uncomment the next line for MS Windows
%     filepath = edfSubfolder + "\" + edfFilename;
%     % uncomment the next line for Linux distributions / Mac OS X
%     % filepath = edfSubfolder + "/" + edfFilename;
% 
%     clear edfSubfolder;
%     clear edfFilename;
%     
%     [hdr, record] = edfread(filepath);
% 
%     wrist_ppg   = record(2,:);
%     f_sample = hdr.frequency(2);
% 
%     time = ((1:size(record, 2)))/f_sample;
% 
%     figure;
%     hold on;
%     subplot(3, 1, 1);
%     plot(time, wrist_ppg);
%     title('Unfiltered');
%     xlabel('time [s]');
%     ylabel('Amplitude [?]');
%     subplot(3, 1, 2);
%     plot(time, PPG_slimBand);
%     title('slimBP filtered');
%     xlabel('time [s]');
%     ylabel('Amplitude [?]');
%     subplot(3, 1, 3);
%     plot(time, PPG_wideBand);
%     title('wideBP filtered');
%     xlabel('time [s]');
%     ylabel('Amplitude [?]');
%     hold off;
%     
%     figure;
%     hold on;
%     plot(time, wrist_ppg);
%     plot(time, PPG_slimBand);
%     plot(time, PPG_wideBand);
%     hold off;

end

function [PP_Temp, AmplitudeCorrectionFactors, PP_PQI, BeatTimes] = ...
    PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, ...
    samplingRate, record)

    % "PPG segmentation, beat localization" in the Papini paper, Figure 1
    
    PP_wideBand = segmentation(PPG_slimBand, PPG_wideBand);
    BeatTimes = beatLocalization(PP_wideBand);
    [PP_Temp, AmplitudeCorrectionFactors, PP_PQI] = pulseNormalization(...
        PP_wideBand);
    
%     % minimum possible time period between heartbeats [s]
%     max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
%     min_RRinterval = 1/max_freq;
%     
%     time = ((1:size(PPG_slimBand, 2)))/samplingRate;
%     time   = time( 1:length(PPG_slimBand) );
%     
%     % find minima
%     PPG_slimBand_inv = PPG_slimBand*-1;
% %     figure;
% %     hold on;
% %     plot(time, PPG_wideBand);
% %     findpeaks(PPG_slimBand, time, 'MinPeakDistance', min_RRinterval);
% %     findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval);
% %     hold off;
%     [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval);
%     BeatTimes = peak_locs;
%         
% %     BeatTimes = % time of detected beats (in seconds after the start of the signal

end

function Temp_Ad = templateCreation(PP_Temp, AmplitudeCorrectionFactors)
    
    % "Template creation" in the Papini paper, Figure 1
    
    % TODO (Ingo)
    
end

function PulseQualityIndex = pulseTemplateComparision(PP_PQI, Temp_Ad)
    
    % "Pulse-Template comparision" in the Papini paper, Figure 1
    
    % TODO (Ingo)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing operations  %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered = bandpassFilter(raw, f_min, f_max, f_sample)
    % "Bandpass filter" in the Papini paper, Figure 1
    f_pass = [f_min f_max] / (0.5 * f_sample);
    [b, a] = butter(3, f_pass, 'bandpass');
    filtered = filter(b, a, raw);
end

function PP = segmentation(PPG_slimBand, PPG_wideBand)
    
    % "Segmentation" in the Papini paper, Figure 1
    
    % TODO (Alex)
    
end

function beatTimes = beatLocalization(PP)
    
    % "Beat Localization" in the Papini paper, Figure 1
    
    % TODO (Alex)
    
end

function [PP_Temp, AmplitudeCorrectionFactors, PP_PQI] = ...
    pulseNormalization(PP)
    
    % "Pulse normalization" in the Papini paper, Figure 1
    
    % TODO (Ingo)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing formulas    %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PP_Temp = PP_Temp(PP)
    % Forumula 9 in the Papini paper
    PP_Temp = ( PP - PPshift(PP) ) / PPamplitude(PP);
end

function PP_PQI = PP_PQI(PP)
    % Formula 10 in the Papini paer
    PP_PQI = PP - PPshift(PP);
end

function PPamplitude = PPamplitude(PP)
    % Formula 7 in the Papini paper
    PPamplitude = abs( max(PP) - min(PP) );
end

function PPshift = PPshift(PP)
    % Formula 8 in the Papini paper
    PPshift = ( max(PP) + min(PP) ) / 2;
end
