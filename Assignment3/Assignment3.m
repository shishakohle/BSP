%% Biosignal Processing - Assingment 2
% Group H: Alexander Neumayer, Ingo Weigel
% 2020-11-08
% Find this code in the project repository on GitHub:
% https://github.com/shishakohle/BSP

clear;
close all;
clc;

%% Choose your ECG .txt file via GUI
[file, path] = uigetfile({'*.txt;'}, 'MultiSelect', 'on');

filepath = path + "/" + file;
ECGdata = importHandscoredRRs(filepath);
clear file path filepath;

% set sample interval and location of first beat in ECG data
sampleInt = 256;
locFirstbeat = 3;

[ECGtimes, ECGintervals] = readECG(ECGdata, sampleInt); % time in seconds after start

function [ECGtimes, ECGintervals] = readECG(ECGdata, sampleInt)

    ECGtimes = cell2mat(ECGdata(3:end, 2)) ./ sampleInt;
    ECGintervals = diff(ECGtimes);

end