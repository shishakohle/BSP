%% Biosignal Processing - Assingment 2
% Group H: Alexander Neumayer, Ingo Weigel
% 2020-11-08
% Find this code in the project repository on GitHub:
% https://github.com/shishakohle/BSP

clear;
close all;
clc;

%% Comment and uncomment according to your ECG .txt file location and your OS
[file, path] = uigetfile({'*.txt;'}, 'MultiSelect', 'on');

filepath = path + "/" + file;
ECGdata = importHandscoredRRs(filepath);
clear file path filepath;

% HandscoredRRs = readECG(ECGdata);

% function [ECGtimes, ECGintervals] = readECG(ECGdata)
% 
% end