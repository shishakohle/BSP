%% REMINDER
% @Alex - don't forget to change "hard coding" of A2 - so that template is
% created for the chosen EDF file

%% Biosignal Processing - Assingment 3
% Group H: Alexander Neumayer, Ingo Weigel
% 2020-11-20
% Find this code in the project repository on GitHub:
% https://github.com/shishakohle/BSP

clear;
close all;
clc;

%% Comment and uncomment according to your EDF file location and your OS

% Define the subfolders for PPG (*.EDF) and ECG (*.TXT) files:
edfSubfolder = "edf";
ecgSubfolder = "ECG";

% Choose on a set of measuerements to be validated:
set = "s1_high_resistance_bike";
% set = "s1_low_resistance_bike";
% set = "s1_walk";
% set = ""; % TODO add further files

% Choose on the path notion of your operating system:

% uncomment the next two lines for MS Windows
% edfFilepath = edfSubfolder + "\" + set + ".edf";
% ecgFilepath = ecgSubfolder + "\" + set + ".txt";

% uncomment the next two lines line for Linux distributions / Mac OS X
edfFilepath = edfSubfolder + "/" + set + ".edf";
ecgFilepath = ecgSubfolder + "/" + set + ".txt";

clear set edfSubfolder ecgSubfolder;

%% Extract raw wrist PPG signal from EDF file

[hdr, rawPPGsignal] = edfread(edfFilepath, 'targetSignals', 'wrist_ppg');
samplingRate        = hdr.samples / hdr.duration;

clear hdr;

%% Extract raw data from ECG file
ECGdata = importHandscoredRRs(ecgFilepath);

% set sample interval and location of first beat in ECG data
sampleInt = 256;
locFirstbeat = 3;

%% Get beat times and intervals from both PPG (according to Assignemnt 2)
%  and ECG raw data

% get PPGtimes and PPGintervals
% For this, execute the essence function of Assignment 2:
[beatTimesAmplitudesPQIs, ~] = PQI(rawPPGsignal, samplingRate);
[PPGbeattimes, PPGbeatintervals] = ...
    readPPG(beatTimesAmplitudesPQIs); % time in seconds after start and intervals in seconds

% get ECGtimes and ECGintervals
[ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat); % time in seconds after start and intervals in seconds

% clear everything we don't need here anymore
clearvars -except ECGbeattimes ECGbeatintervals PPGbeattimes PPGbeatintervals rawPPGsignal PPGsignal ECGdata samplingRate;

%% Call validate function

% real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes);

% doesnt work because of different dimensions - for loop needed
% function real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes)
    
    % pulse transit time calculation/estimation
    avgPWV = 6.84; % average pulse wave velocity [m/s] of healthy persons Diaz et. al. "Reference Values of Pulse Wave Velocity in Healthy People from an Urban and Rural Argentinean Population", International Journal of Hypertension, vol. 2014, Article ID 653239, 7 pages, 2014. https://doi.org/10.1155/2014/653239
    avgHeightFM = 1.66; % average body height of west european women [m] https://www.worlddata.info/average-bodyheight.php
    avgHeightM = 1.8; % average body height of west european men [m] https://www.worlddata.info/average-bodyheight.php
    avgHeight = (avgHeightFM + avgHeightM) / 2;
    avgLengthHearttoFinger = avgHeight / 2; % average length from the heart to the fingertip of west europeans [m] (https://www.scientificamerican.com/article/human-body-ratios/)
    avgPTT = avgLengthHearttoFinger / avgPWV; % average PTT from heart to fingertip in [s]
    
    % ECG beats corrected to PPG signal by adding avgPTT
    ECGplusPTTtimes = ECGbeattimes + avgPTT;
    Tolerance = 0.1; % pecentage
    timeTolerance = avgPTT + avgPTT * Tolerance;

    [TP, locTPandFN] = ismembertol(PPGbeattimes, ECGplusPTTtimes, timeTolerance, 'DataScale', 1); % TP = logical vector containing 1 at the index for all true positives and 0 for all false positives / FN = vector containing index of the elements where a beat in PPG was found in ECG and 0 where it is a false negative
    
    % calculate specificity
    nTP = sum(TP);
    nFP = size(PPGbeattimes, 1) - nTP;
    nFN = size(ECGbeattimes, 1) - size(find(~locTPandFN), 1); 
    Specificity = nTP / (nTP + nFN);
    
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

% end

%% Function declaration
function [ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat)

    ECGbeattimes = cell2mat(ECGdata(locFirstbeat:end, 2)) ./ sampleInt;
    ECGbeatintervals = diff(ECGbeattimes);

end

function [PPGbeattimes, PPGbeatintervals] = readPPG(PPGdata)

    PPGbeattimes = transpose(cell2mat(PPGdata(1, 1)));
    PPGbeatintervals = diff(PPGbeattimes);

end