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

% choose on an EDF subfolder:
% edfSubfolder = ".";
edfSubfolder = "edf";

% choose on an EDF file:
edfFilename  = "s1_high_resistance_bike.edf";
% edfFilename  = ""; % TODO add further files
% edfFilename  = "";

% choose on the path notion of your operating system:
% uncomment the next line for MS Windows
% filepath = edfSubfolder + "\" + edfFilename;
% uncomment the next line for Linux distributions / Mac OS X
filepath = edfSubfolder + "/" + edfFilename;

clear edfSubfolder;
clear edfFilename;

%% Extract raw wrist PPG signal from EDF file

[hdr, rawPPGsignal] = edfread(filepath, 'targetSignals', 'wrist_ppg');
samplingRate        = hdr.samples / hdr.duration;

%% Execute the essence function of Assignment 2

[beatTimesAmplitudesPQIs, ~] = PQI(rawPPGsignal, ...
    samplingRate);

%% Choose your ECG .txt file via GUI
[file, path] = uigetfile({'*.txt;'}, 'MultiSelect', 'on');

filepath = path + "/" + file;
ECGdata = importHandscoredRRs(filepath);

% set sample interval and location of first beat in ECG data
sampleInt = 256;
locFirstbeat = 3;

% get ECGtimes and ECGintervals
[ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat); % time in seconds after start and intervals in seconds

% get PPGtimes and PPGintervals
[PPGbeattimes, PPGbeatintervals, PPGsignal] = readPPG(beatTimesAmplitudesPQIs); % time in seconds after start and intervals in seconds

% clear everything we don't need here anymore
clearvars -except ECGbeattimes ECGbeatintervals PPGbeattimes PPGbeatintervals rawPPGsignal PPGsignal ECGdata samplingRate;

%% Call validate function

% real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes);

% doesnt work because of different dimensions - for loop needed
% function real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes)
    
    % pulse transit time calculation/estimation
    avgPWV = 6.84; % average pulse wave velocity [m/s] of healthy persons Díaz et. al. "Reference Values of Pulse Wave Velocity in Healthy People from an Urban and Rural Argentinean Population", International Journal of Hypertension, vol. 2014, Article ID 653239, 7 pages, 2014. https://doi.org/10.1155/2014/653239
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
    time_PPG = ((1:size(PPGsignal, 1)))/samplingRate;

    ECGvect = zeros(size(ECGbeattimes, 1), 1);

    figure;
    hold on;
    findpeaks(PPGsignal, time_PPG, 'MinPeakDistance', min_RRinterval);
    plot(ECGbeattimes, ECGvect, 'r*');
    plot(ECGplusPTTtimes, ECGvect, 'g*');
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