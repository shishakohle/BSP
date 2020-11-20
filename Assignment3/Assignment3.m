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
[beatTimesAmplitudesPQIs, ~]     = PQI(rawPPGsignal, samplingRate);
[PPGbeattimes, PPGbeatintervals] = readPPG(beatTimesAmplitudesPQIs); % time in seconds after start and intervals in seconds

% get ECGtimes and ECGintervals
[ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat); % time in seconds after start and intervals in seconds

% clear everything we don't need here anymore
clearvars -except ECGbeattimes ECGbeatintervals PPGbeattimes PPGbeatintervals rawPPGsignal PPGsignal ECGdata samplingRate;

%% Call validate function

% real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes);

% doesnt work because of different dimensions - for loop needed
% function real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes)
    
    % pulse transit time calculation/estimation
    avgPWV = 6.84; % average pulse wave velocity [m/s] of healthy persons D�az et. al. "Reference Values of Pulse Wave Velocity in Healthy People from an Urban and Rural Argentinean Population", International Journal of Hypertension, vol. 2014, Article ID 653239, 7 pages, 2014. https://doi.org/10.1155/2014/653239
    avgHeightFM = 1.66; % average body height of west european women [m]
    avgHeightM = 1.8; % average body height of west european men [m]
    avgHeight = (avgHeightFM + avgHeightM) / 2;
    avgLengthHearttoFinger = avgHeight / 2; % average length from the heart to the fingertip of west europeans [m] (https://www.scientificamerican.com/article/human-body-ratios/)
    avgPTT = avgLengthHearttoFinger / avgPWV; % average PTT from heart to fingertip in [s]
    
    ECGplusPTTtimes = ECGbeattimes + avgPTT;
    timeTolerance = 0.5; % pecentage
    timeLowerLimit = ECGplusPTTtimes - (avgPTT + avgPTT * timeTolerance);
    timeUpperLimit = ECGplusPTTtimes + (avgPTT + avgPTT * timeTolerance);
    count = 1;
    
    % Plot PPG and ECG

    max_freq = 210/60; % maximum expectable frequency is 210 bpm -> converted to Hz
    min_RRinterval = 1/max_freq;
%     time_PPG = ((1:size(PPGsignal, 1)))/samplingRate;

    ECGvect = zeros(size(ECGbeattimes, 1), 1);
    PPGvect = zeros(size(PPGbeattimes, 1), 1);

    figure;
    hold on;
%     findpeaks(PPGsignal, time_PPG, 'MinPeakDistance', min_RRinterval);
    plot(ECGbeattimes, ECGvect, 'r*');
    plot(ECGplusPTTtimes, ECGvect, 'g*');
    plot(PPGbeattimes, PPGvect, 'b*');
    hold off;

    for i = 1 : size(PPGbeattimes, 1)
        if PPGbeattimes(i) >= timeLowerLimit(i) && PPGbeattimes(i) <= timeUpperLimit(i)
        real_beats(count) = PPGbeattimes(i);
        count = count + 1;
        end
    end

% end

%% Function declaration
function [ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat)

    ECGbeattimes = cell2mat(ECGdata(locFirstbeat:end, 2)) ./ sampleInt;
    ECGbeatintervals = diff(ECGbeattimes);

end

function [PPGbeattimes, PPGbeatintervals, PPGsignal] = readPPG(PPGdata)

    PPGbeattimes = transpose(cell2mat(PPGdata(1, 1)));
    PPGbeatintervals = diff(PPGbeattimes);
    PPGsignal = transpose(cell2mat(PPGdata(4, 1)));

end