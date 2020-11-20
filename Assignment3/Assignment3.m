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
clearvars -except ECGbeattimes ECGbeatintervals PPGbeattimes PPGbeatintervals;

%% Call validate function

Sensitvity = Validate(PPGbeattimes, PPGbeatintervals, ECGbeattimes, ECGbeatintervals);

%% Function declaration
function [ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat)

    ECGbeattimes = cell2mat(ECGdata(locFirstbeat:end, 2)) ./ sampleInt;
    ECGbeatintervals = diff(ECGbeattimes);

end

function [PPGbeattimes, PPGbeatintervals] = readPPG(PPGdata)

    PPGbeattimes = transpose(cell2mat(PPGdata(1, 1)));
    PPGbeatintervals = diff(PPGbeattimes);

end