%% Function PQI() as a solution to Assignment 3

function [Sensitivity, ECGbeatintervalsforAnalysis, PPGbeatintervalsforAnalysis, nTP, nFP, nFN, beatplotFigure, scatterplotFigure] ...
    = Validate(PPGbeattimes, PPGbeatintervals, ECGbeattimes, ...
    ECGbeatintervals)
    
    % pulse transit time calculation/estimation
    avgPWV = 6.84; % average pulse wave velocity [m/s] of healthy persons Diaz et. al. "Reference Values of Pulse Wave Velocity in Healthy People from an Urban and Rural Argentinean Population", International Journal of Hypertension, vol. 2014, Article ID 653239, 7 pages, 2014. https://doi.org/10.1155/2014/653239
    avgHeightFM = 1.66; % average body height of west european women [m] https://www.worlddata.info/average-bodyheight.php
    avgHeightM = 1.8; % average body height of west european men [m] https://www.worlddata.info/average-bodyheight.php
    avgHeight = (avgHeightFM + avgHeightM) / 2;
    avgLengthHearttoFinger = avgHeight / 2; % average length from the heart to the fingertip of west europeans [m] (https://www.scientificamerican.com/article/human-body-ratios/)
    avgPTT = avgLengthHearttoFinger / avgPWV; % average PTT from heart to fingertip in [s]
    
    % ECG beats corrected to PPG signal by adding avgPTT
    ECGplusPTTtimes = ECGbeattimes + avgPTT;
    Tolerance = 0.15; % 15 pecent tolerance - determined empirically
    timeTolerance = avgPTT + avgPTT * Tolerance;
    timeTolerance = timeTolerance + 0.125 / 2; % account for the pulse transit time variations in PPG according to Foo et al 2005

    [TP, locTPandFN] = ismembertol(PPGbeattimes, ECGplusPTTtimes, timeTolerance, 'DataScale', 1); % TP = logical vector containing 1 at the index for all true positives and 0 for all false positives / locTPandFN = vector containing index of the elements where a beat in PPG was found in ECG and 0 where it is a false negative
    
    % calculate sensitivity of beat detection
    nTP = sum(TP);
    nFP = size(PPGbeattimes, 1) - nTP;
    nFN = size(ECGbeattimes, 1) - size(find(locTPandFN), 1); 
    Sensitivity = nTP / (nTP + nFN);
    PPGtoECGbeat_ratio = length(PPGbeattimes) / length(ECGbeattimes);
    
    % make plot for graphical investigation - green ECG beats with PPT
    % shift / blue PPG beats / black TP PPG beats squares are on 1 - so if
    % green star is not in a balck square its a FP
    ECGplusPTTtimes_lowerlim = ECGplusPTTtimes - timeTolerance;
    ECGplusPTTtimes_upperlim = ECGplusPTTtimes + timeTolerance;
    
    ECGvect = zeros(size(ECGbeattimes, 1), 1);
    PPGvect = zeros(size(PPGbeattimes, 1), 1);
    beatplotFigure = figure;
    beatplotFigure.Visible='on';
    hold on;
    plot(ECGbeattimes, ECGvect, 'r*');
    plot(ECGplusPTTtimes, ECGvect, 'g*');
    plot(ECGplusPTTtimes_lowerlim, ECGvect, 'g>');
    plot(ECGplusPTTtimes_upperlim, ECGvect, 'g<');
    plot(PPGbeattimes, PPGvect, 'b*');
    plot(PPGbeattimes, TP, 'ks');
    ylim([0 1.5]);
    legend('ECG beats','Adjusted ECG beats', 'time tolerance lower limit', 'time tolerance upper limit', 'PPG beats', 'false positive PPG beats'); 
    hold off;

    % make bland-altman and correlation plots for inter beat intervals of
    % the TP beats
    locTPinECGbeats = locTPandFN;
    locTPinECGbeats(~locTPandFN) = [];
    if locTPinECGbeats(end) > length(ECGbeatintervals)
        
        locTPinECGbeats(end) = [];
        
    end
    ECGbeatintervalsforAnalysis = ECGbeatintervals(locTPinECGbeats); % could lead to problem if index in locTPinECGbeats > size of ECGbeatintervals
    if length(TP) > length(PPGbeatintervals)
        
        TP(end) = [];
        
    end
    PPGbeatintervalsforAnalysis = PPGbeatintervals(TP);
    
    if length(PPGbeatintervalsforAnalysis) > length(ECGbeatintervalsforAnalysis)
        
        PPGbeatintervalsforAnalysis(end) = [];
        
    end
    
    % BA plot paramters
    tit = 'BeatInterval analysis'; % figure title
    label = {'ECG beat intervals', 'PPG beat intervals', 'seconds'}; % Names of data sets

    [~, scatterplotFigure, ~] = BlandAltman(ECGbeatintervalsforAnalysis, PPGbeatintervalsforAnalysis, label, tit);
    scatterplotFigure.Visible='off';
end