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
clear file path filepath;

% set sample interval and location of first beat in ECG data
sampleInt = 256;
locFirstbeat = 3;

% get ECGtimes and ECGintervals
[ECGtimes, ECGintervals] = readECG(ECGdata, sampleInt, locFirstbeat); % time in seconds after start and intervals in seconds

% get PPGtimes and PPGintervals
[PPGtimes, PPGintervals] = readPPG(beatTimesAmplitudesPQIs); % time in seconds after start and intervals in seconds

% clear everything we don't need here anymore
clear ECGdata hdr locFirstbeat rawPPGsignal sampleInt samplingRate;

%% Call validate function

%% Function declaration
function [ECGtimes, ECGintervals] = readECG(ECGdata, sampleInt, locFirstbeat)

    ECGtimes = cell2mat(ECGdata(locFirstbeat:end, 2)) ./ sampleInt;
    ECGintervals = diff(ECGtimes);

end

function [PPGtimes, PPGintervals] = readPPG(PPGdata)

    PPGtimes = transpose(PPGdata(1, :));
    PPGintervals = diff(PPGtimes);

end