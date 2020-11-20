%% Function PQI() as a solution to Assignment 3

function Sensitvity = Validate(PPGbeattimes, PPGbeatintervals, ECGbeattimes, ECGbeatintervals)
    
    % pulse transit time calculation/estimation
    avgPWV = 6.84; % average pulse wave velocity [m/s] of healthy persons Diaz et. al. "Reference Values of Pulse Wave Velocity in Healthy People from an Urban and Rural Argentinean Population", International Journal of Hypertension, vol. 2014, Article ID 653239, 7 pages, 2014. https://doi.org/10.1155/2014/653239
    avgHeightFM = 1.66; % average body height of west european women [m] https://www.worlddata.info/average-bodyheight.php
    avgHeightM = 1.8; % average body height of west european men [m] https://www.worlddata.info/average-bodyheight.php
    avgHeight = (avgHeightFM + avgHeightM) / 2;
    avgLengthHearttoFinger = avgHeight / 2; % average length from the heart to the fingertip of west europeans [m] (https://www.scientificamerican.com/article/human-body-ratios/)
    avgPTT = avgLengthHearttoFinger / avgPWV; % average PTT from heart to fingertip in [s]
    
    % ECG beats corrected to PPG signal by adding avgPTT
    ECGplusPTTtimes = ECGbeattimes + avgPTT;
    Tolerance = 0.3; % 30 pecent tolerance - determined empirically
    timeTolerance = avgPTT + avgPTT * Tolerance;

    [TP, locTPandFN] = ismembertol(PPGbeattimes, ECGplusPTTtimes, timeTolerance, 'DataScale', 1); % TP = logical vector containing 1 at the index for all true positives and 0 for all false positives / locTPandFN = vector containing index of the elements where a beat in PPG was found in ECG and 0 where it is a false negative
    numberuniqueentriesinlocTPandFN = length(unique(locTPandFN));
    
    % tried an approach to find a tolerance value which fits so that only unique PPG
    % beats are found in ECG data - but it doesnt work. I'm too tired at this
    % point
%     while numberuniqueentriesinlocTPandFN < length(PPGbeattimes)
%         
%         Tolerance = Tolerance - 0.1;
%         timeTolerance = avgPTT + avgPTT * Tolerance;
%         [TP, locTPandFN] = ismembertol(PPGbeattimes, ECGplusPTTtimes, timeTolerance, 'DataScale', 1); % TP = logical vector containing 1 at the index for all true positives and 0 for all false positives / locTPandFN = vector containing index of the elements where a beat in PPG was found in ECG and 0 where it is a false negative
%         numberuniqueentriesinlocTPandFN = length(unique(locTPandFN));
% 
%     end
    
    % calculate sensitivity of beat detection
    nTP = sum(TP);
    nFP = size(PPGbeattimes, 1) - nTP;
    nFN = size(ECGbeattimes, 1) - size(find(locTPandFN), 1); 
    Sensitvity = nTP / (nTP + nFN);
    
    % make plot for graphical investigation - green ECG beats with PPT
    % shift / blue PPG beats / black TP PPG beats squares are on 1 - so if
    % green star is not in a balck square its a FP
    ECGvect = zeros(size(ECGbeattimes, 1), 1);
    PPGvect = zeros(size(PPGbeattimes, 1), 1);
    figure;
    hold on;
%     plot(ECGbeattimes, ECGvect, 'r*');
    plot(ECGplusPTTtimes, ECGvect, 'g*');
    plot(PPGbeattimes, PPGvect, 'b*');
    plot(PPGbeattimes, TP, 'ks');
    ylim([0 1.5]);
    hold off;

    % make bland-altman and correlation plots for inter beat intervals of
    % the TP beats
    locTPinECGbeats = locTPandFN;
    locTPinECGbeats(~locTPandFN) = [];
    ECGbeatintervalsforAnalysis = ECGbeatintervals(locTPinECGbeats); % could lead to problem if index in locTPinECGbeats > size of ECGbeatintervals
    PPGbeatintervalsforAnalysis = PPGbeatintervals(TP);
    
    % BA plot paramters
    corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
    BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot

    BlandAltman(ECGbeatintervalsforAnalysis, PPGbeatintervalsforAnalysis, 'BeatIntervals', 'BeatInterval analysis');
end