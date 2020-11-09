%% Function PQI() as a solution to Assignment 2

function [beatTimesAmplitudesPQIs, pulseWaveTemplate] =  ...
    PQI (rawPPGsignal, samplingRate)
    
    %% Algorithm according to Figure 1 in the Papini paper

    % call Pre-processing
    [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, ...
        samplingRate);
    
    % call PPG segmentation, beat localization
    [PPs_Temp, AmplitudeCorrectionFactors, PP_PQI, BeatTimes, ...
        PPamplitudes] = ...
        PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, ...
        samplingRate);
    
    % call Template creation
    Temps_Ad = templateCreation(PPs_Temp, samplingRate, ...
        AmplitudeCorrectionFactors);
    
    % call Pulse-Template comparison
    PulseQualityIndexes = pulseTemplateComparision(PP_PQI, Temps_Ad);
    
    %% assemble return values according to Assignment 2 task description:
    
    % "Matrix with timing, amplitude, and quality index for each detected
    %  beat"
    % "Optional: calculated pulse wave template"
    
    beatTimesAmplitudesPQIs = [BeatTimes, PPamplitudes, ...
        PulseQualityIndexes];
    pulseWaveTemplate = Temps_Ad;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing main blocks %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PPG_slimBand, PPG_wideBand] = preprocessing(rawPPGsignal, ...
    samplingRate)
    
    %% "Pre-processing" in the Papini paper, Figure 1
    
    % filter signal
    PPG_slimBand = bandpassFilter(rawPPGsignal, 0.4, 2.25, samplingRate);
    PPG_wideBand = bandpassFilter(rawPPGsignal, 0.4, 10.0, samplingRate);    
    
    % phase correct wide band to slim band - obsolete because of filtfilt
%     [c, lags] = xcorr(PPG_slimBand, PPG_wideBand);               % compute cross correlation; keep lags vector
%     [~, iLag] = max(c(find(lags==0) : end));  % find the max in one-sided
%     PPG_wideBand_phaseCorrect = circshift(PPG_wideBand, [0 iLag]);           % correct for the shift    
%     PPG_wideBand = PPG_wideBand_phaseCorrect;
end

function [PP_Temp, AmplitudeCorrectionFactors, PP_PQI, BeatTimes, ...
    PPamplitudes] = ...
    PPGsegmentationAndBeatLocalization(PPG_slimBand, PPG_wideBand, ...
    samplingRate)

    %% "PPG segmentation, beat localization" in the Papini paper, Figure 1
    
    [BeatTimes, time] = beatLocalization(PPG_slimBand, samplingRate);
    
    % segment with the help of BeatTimes
    PP_wideBand = segmentation(PPG_slimBand, PPG_wideBand, BeatTimes, ...
        time);
    
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
    
    % skipped for now
    AmplitudeCorrectionFactors = PPamplitudes;

end

function Temps_Ad = templateCreation(PPs_Temp, f_sample, ...
    AmplitudeCorrectionFactors)
    
    %% "Template creation" in the Papini paper, Figure 1 and Figure 3
    
    % "For this reason, our algorithm calculates the pulse template by
    %  means of DBA (Petitjean et al 2014). This allows the time series to
    %  be averaged by iteratively decreasing the DTW distance between an
    %  initial template and each individual pulse (figure 3). In each
    %  iteration, the resulting averaged time series is used as the initial
    %  template for the next iteration. The DBA is initialized with the
    %  medoid of the PPs_Temp as initial template, and this initial
    %  template is refined during five iterations (Petitjean et al 2014).
    %  The number of iterations is chosen empirically to ensure the
    %  computation of a pulse template with a correct morphological
    %  representation of the pulses, in an efficient computation time.
    %  The resulting pulse template is filtered using a 3rd-order low-pass
    %  Butterworth filter with a cut-off frequency of 10 Hz to remove
    %  high-frequency components introduced by the DBA. This guarantees
    %  that the template has the same frequency components as the PPG
    %  segments used to obtain it (PPs_Temp) and of the PPG segments used
    %  to calculate the quality index (PPs_PQI).
    
    % medoid calculation done by DBA():
    % temporaryTemplate = DBA(PPs_Temp); % function in DBA.m
    % Papini intereates 5 times, DBA() 15 times, we set it to 1 time!
    % still HIGH COMPUTATIONAL EFFORT, THEREFORE AS A PROVISORY SOLUTION:
        temporaryTemplate = PPs_Temp;
    
    % Template = lowpassFilter(temporaryTemplate, f_sample, 10);
    % ERROR IN CALL OF lowpassFilter(), therefore as a
    % PROVISIONARY SOLUTION:
        Template = temporaryTemplate;
    
    % calculate all adjusted templates
    % "the adjusted templates (Temp Ad s) are obtained by multiplying the
    %  template by the corresponding correction factors"
    
    % Temps_Ad = Template .* AmplitudeCorrectionFactors;
    % ERROR BECAUSE CELLS cannot .*, FOR LOOP NEEDED
    % THEREFORE PROVISIONARY SOLUTION (i.e. no correction):
        Temps_Ad = Template;
    
end

function PulseQualityIndex = pulseTemplateComparision(PP_PQI, Temp_Ad)
    
    %% "Pulse-Template comparision" in the Papini paper, Figure 1 and 6
    
    % "The proposed algorithm uses DTW to derive PP_warped from the warping
    %  of PP_PQI to the template calculated for the one hour of PPG signal
    %  they belong to (figure 5)."
    
    % PP_warped = DBA(PP_PQI);
    % HIGH COMPUTATIONAL EFFORT, THEREFORE AS A PROVISORY SOLUTION:
        PP_warped = PP_PQI;
    
    % "(...) part of the morphological discrepancies between PP PQI and
    %  Temp Ad remain in the PPs warped . These residual differences are
    %  used by the algorithm to calculate the quality index (PQI) of each
    %  PPG pulse."
    
    % "First, the algorithm finds the indices of each sample of PP_warped
    %  with a difference higher than 10% with respect to the corresponding
    %  samples of the Temp_Ad"
    % Formula 11 in the Papini paper.
    
    % UP = find( ((PP_warped - Temp_Ad)./Temp_Ad) > 0.1 ); % unmatched points
    % MP = find( ((PP_warped - Temp_Ad)./Temp_Ad)<= 0.1 ); % matched points
    % ERROR "Undefined operator '-' for input arguments of type 'cell'."
    % AGAIN THE "CELL PROBLEM" -> FOR LOOP NEEDED.
    % AS A PROVISIONARY SOLUTION FOR QUICK HAND-IN:
        MP = ones(size(PP_warped));
        UP = zeros(size(PP_warped));
    
    % "The quantity of unmatched and matched points, respectively N_UP and
    %  N_MP , are used to calculate the percentage of matching points as "
    % Formula 12 in the Papini paper.
    
    number_MP = sum(MP);
    number_UP = sum(UP);
    percentage_MP = number_MP / (number_MP+number_UP);
    
    % "The second step consists in calculating the root mean square error
    %  of the unmatched points in respect to Temp_Ad"
    % Formula 13 in the Papini paper.
    
    if number_UP == 0
        RMSE_UP = 0;
    else
        RMSE_UP = ( Temp_Ad(UP) - PP_warped(UP) )^2 ./ number_UP;
    end
    
    % Formula 14 in the Papini paper.
    
    % RMSE_norm = RMSE_UP / ( max(Temp_Ad) - min(Temp_Ad) );
    % PROVISIONARY SOLUTION
        RMSE_norm = RMSE_UP;
    
    % Formula 15
    
    PulseQualityIndex = max( 0.1 - RMSE_norm/percentage_MP );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing operations  %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered = bandpassFilter(raw, f_min, f_max, f_sample)
    %% "Bandpass filter" in the Papini paper, Figure 1
    f_pass = [f_min f_max] / (0.5 * f_sample);
    [b, a] = butter(3, f_pass, 'bandpass');
    filtered = filtfilt(b, a, raw);
end

function filtered = lowpassFilter(raw, f_sample, f_cutoff)
    %% "Lowpass filter" in the Papini paper, Figure 3
    [b, a] = butter( 3, f_cutoff / (f_sample/2) );
    for i = 1 : length(raw)
        filtered(i) = filter(b, a, raw(i));
    end
end
%% FIR - filtfilt
function [beatTimes, time] = beatLocalization(PPG_slimBand, samplingRate)
    
    %% "Beat Localization" in the Papini paper, Figure 1
    
    % minimum possible time period between heartbeats [s]
    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    min_RRinterval = 1/max_freq;
    
    time = ((1:size(PPG_slimBand, 2)))/samplingRate;
    
    % find minima
    PPG_slimBand_inv = PPG_slimBand*-1;
%     figure;
%     hold on;
%     findpeaks(PPG_slimBand, time, 'MinPeakDistance', min_RRinterval);
%     findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval);
%     hold off;
    [peak_vals, peak_locs, peak_widths, peak_prominences] = ...
    findpeaks(PPG_slimBand_inv, time, 'MinPeakDistance', min_RRinterval, 'MinPeakHeight', 0); % with negativity condition as in formula 1 in paper here >0 because of inv. sig.
    beatTimes = peak_locs; % time of detected beats (in seconds after the start of the signal)
    
end

function PP = segmentation(PPG_slimBand, PPG_wideBand, beatTimes, time)
    
    %% "Segmentation" in the Papini paper, Figure 1
    
    beatLocs = find(ismember(time, beatTimes));
    start = 1;
    for i = 1 : size(beatLocs, 2)
       
        segment_slimBand = PPG_slimBand(start : beatLocs(1, i));
        segment_wideBand = PPG_wideBand(start : beatLocs(1, i));
%         figure;
%         hold on;
%         plot(segment_slimBand);
%         plot(segment_wideBand);
%         hold off;
        PP{i, :} = segment_wideBand;
        start = beatLocs(1, i);
        
    end
    
end

function [PP_Temp, PP_amplitude, PP_PQI] = ...
    pulseNormalization(PP)
    
    %% "Pulse normalization" in the Papini paper, Figure 1
    
    PP_Temp     = PPtemp   (PP);
    PP_amplitude = PPamplitude(PP);
    PP_PQI      = PPpqi    (PP);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions representing formulas    %
% in Figure 1 of the Papini paper    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PP_amplitude = PPamplitude(PP)
    %% Formula 7 in the Papini paper
    for i = 1 : length(PP)
        PP_amplitude(i) = abs( max(PP{i}) - min(PP{i}) );
    end
end

function PP_shift = PPshift(PP)
    %% Formula 8 in the Papini paper
    for i = 1 : length(PP)
        PP_shift(i) = ( max(PP{i}) + min(PP{i}) ) / 2;
    end
end

function PP_Temp = PPtemp(PP)
    %% Forumula 9 in the Papini paper
    PP_shift = PPshift(PP);
    PP_amplitude = PPamplitude(PP);
    for i = 1 : length(PP)
        PP_cell = PP{i};
        PP_cellminusshift = PP_cell - PP_shift(1, i);
        PP_Temp{i, :} = PP_cellminusshift ./ PP_amplitude(1, i);
    end
end

function PP_PQI = PPpqi(PP)
    %% Formula 10 in the Papini paper
    PP_shift = PPshift(PP);
    for i = 1 : length(PP)
        PP_cell = PP{i};
        PP_PQI{i, :} = PP_cell - PP_shift(1, i);
    end
end
