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
    
    % filter signal
    PPG_slimBand = bandpassFilter(rawPPGsignal, 0.4, 2.25, samplingRate);
    PPG_wideBand = bandpassFilter(rawPPGsignal, 0.4, 10.0, samplingRate);
    
    % phase correct wide band to slim band
    [c, lags] = xcorr(PPG_slimBand, PPG_wideBand);               % compute cross correlation; keep lags vector
    [~, iLag] = max(c(find(lags==0) : end));  % find the max in one-sided
    PPG_wideBand_phaseCorrect = circshift(PPG_wideBand, [0 iLag]);           % correct for the shift
    
    
    % optionally plot each signal - only for code testing, exclude in the
    % final version
    
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
    plot(time, wrist_ppg);
    plot(time, PPG_slimBand);
    plot(time, PPG_wideBand);
    plot(time, PPG_wideBand_phaseCorrect, 'g');
    hold off;
    
    PPG_wideBand = PPG_wideBand_phaseCorrect;

end

function [PP_Temp, AmplitudeCorrectionFactors, PP_PQI, BeatTimes] = ...
    PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, ...
    samplingRate, record)

    % "PPG segmentation, beat localization" in the Papini paper, Figure 1
    
    [BeatTimes, time] = beatLocalization(PPG_slimBand, samplingRate);    
    PP_wideBand = segmentation(PPG_slimBand, PPG_wideBand, BeatTimes, time); % �ber beat times segementieren
    [PP_Temp, PPamplitudes, PP_PQI] = pulseNormalization(...
        PP_wideBand);
    
    % "To derive the correction factors, a time series comprising all the
    %  amplitudes obtained in (7) is stored. First, the algorithm removes
    %  from the amplitude time series all elements that have a value 50%
    %  higher or lower than the previous or the following value. Then the
    %  clipped amplitude time series is interpolated at 4 Hz using a cubic
    %  spline interpolation and filtered with a 3rd-order zero-phase
    %  low-pass Butterworth filter with a cut-off frequency of 1.5 Hz.
    %  (...) Finally, the filtered signal is resampled at the same time
    %  locations of the original amplitude time series."

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

function [beatTimes, time] = beatLocalization(PPG_slimBand, samplingRate)
    
    % "Beat Localization" in the Papini paper, Figure 1
    
    % minimum possible time period between heartbeats [s]
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    min_RRinterval = 1/max_freq;
    
    time = ((1:size(PPG_slimBand, 2)))/samplingRate;
    
    % find minima
    PPG_slimBand_inv = PPG_slimBand*-1;
%     figure;
%     hold on;
%     plot(time, PPG_wideBand);
%     findpeaks(PPG_slimBand, time, 'MinPeakDistance', min_RRinterval);
%     findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval);
%     hold off;
    [peak_vals, peak_locs, peak_widths, peak_prominences] = findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval);
    beatTimes = peak_locs; % time of detected beats (in seconds after the start of the signal)
    
end

function PP = segmentation(PPG_slimBand, PPG_wideBand, beatTimes, time)
    
    % "Segmentation" in the Papini paper, Figure 1
    
    beatLocs = find(ismember(time, beatTimes));
    start = 1;
    for i = 1 : size(beatLocs, 2)
       
        segment_slimBand = PPG_slimBand(start : beatLocs(1, i));
        segment_wideBand = PPG_wideBand(start : beatLocs(1, i));
        figure;
        hold on;
        plot(segment_slimBand);
        plot(segment_wideBand);
        hold off;
        PP{i, :} = segment_wideBand;
        start = beatLocs(1, i);
        
    end
    
end

function [PP_Temp, PP_Amplitude, PP_PQI] = ...
    pulseNormalization(PP)
    
    % "Pulse normalization" in the Papini paper, Figure 1
    
    PP_Temp   = PP_Temp(PP);
    PPamplitude = PPamplitude(PP);
    PP_PQI    = PP_PQI (PP);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing formulas    %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PPamplitude = PPamplitude(PP)
    % Formula 7 in the Papini paper
    PPamplitude = abs( max(PP) - min(PP) );
end

function PPshift = PPshift(PP)
    % Formula 8 in the Papini paper
    PPshift = ( max(PP) + min(PP) ) / 2;
end

function PP_Temp = PP_Temp(PP)
    % Forumula 9 in the Papini paper
    PP_Temp = ( PP - PPshift(PP) ) / PPamplitude(PP);
end

function PP_PQI = PP_PQI(PP)
    % Formula 10 in the Papini paer
    PP_PQI = PP - PPshift(PP);
end
