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
Datasets = ["s1_low_resistance_bike", "s1_walk", "s2_high_resistance_bike", ... 
    "s2_low_resistance_bike", "s2_walk", "s3_high_resistance_bike", ...
    "s3_low_resistance_bike", "s3_run", "s3_walk", "s4_run", ...
    "s5_low_resistance_bike", "s5_run", "s6_low_resistance_bike", ...
    "s6_run", "s6_walk", "s8_run", "s8_walk", "s9_walk"];

% Choose on the path notion of your operating system:

% uncomment the next two lines for MS Windows
% edfFilepath = edfSubfolder + "\" + set + ".edf";
% ecgFilepath = ecgSubfolder + "\" + set + ".txt";

% uncomment the next two lines line for Linux distributions / Mac OS X
edfFilepath = edfSubfolder + "/" + Datasets + ".edf";
ecgFilepath = ecgSubfolder + "/" + Datasets + ".txt";

clear set edfSubfolder ecgSubfolder;

%% Get Beattimes and Intervals

for i = 1 : length(Datasets)
   
    [PPGbeattimes{i}, PPGbeatintervals{i}, ECGbeattimes{i}, ECGbeatintervals{i}] = ...
    getBeattimesandIntervals(edfFilepath(i), ecgFilepath(i));
    
end

%% Call validate function and calculate mean sensitivity and countinous output metrics

All_ECGbeatintervalsforAnalysis = [];
All_PPGbeatintervalsforAnalysis = [];
for i = 1 : length(Datasets)

    [Sensitivity{i}, ECGbeatintervalsforAnalysis{i}, PPGbeatintervalsforAnalysis{i}, nTP{i}, nFP{i}, nFN{i}, beatplotFigure, scatterplotFigure] = ...
        Validate(PPGbeattimes{i}, PPGbeatintervals{i}, ECGbeattimes{i}, ...
        ECGbeatintervals{i});

%     % print results to console:
%     display ("Sensitvity: "+ Sensitvity);
%     display ("True Positives: " + nTP);
%     display ("False Positives: " + nFP);
%     display ("False Negatives: " + nFN);
% 
%     % you may want to show one of the returned plots, e.g. the scatterplot:
    scatterplotFigure.Visible='off';
    
    start = length(All_ECGbeatintervalsforAnalysis)+1;
    ending = start + length(ECGbeatintervalsforAnalysis{i})-1;
    All_ECGbeatintervalsforAnalysis(start : ending) = ECGbeatintervalsforAnalysis{i};
    
    start = length(All_PPGbeatintervalsforAnalysis)+1;
    ending = start + length(PPGbeatintervalsforAnalysis{i})-1;
    All_PPGbeatintervalsforAnalysis(start : ending) = PPGbeatintervalsforAnalysis{i};
    
end

% calculate mean sensitivity
MeanSensitivity = mean(cell2mat(Sensitivity));

% BA plot paramters
tit = 'BeatInterval analysis'; % figure title
gnames = {ECGbeatintervals, PPGbeatintervals}; % names of groups in data {dimension 1 and 2}
label = {'ECG beat intervals', 'PPG beat intervals', 'seconds'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks', 'CV'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
colors = 'br';      % character codes


% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman(All_ECGbeatintervalsforAnalysis, All_PPGbeatintervalsforAnalysis,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'symbols','Num', 'showFitCI',' on');


%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions %%
%%%%%%%%%%%%%%%%%%%%%%

function [PPGbeattimes, PPGbeatintervals, ECGbeattimes, ECGbeatintervals] = getBeattimesandIntervals(edfFilepath, ecgFilepath)

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

end

function [ECGbeattimes, ECGbeatintervals] = readECG(ECGdata, sampleInt, locFirstbeat)

    ECGbeattimes = cell2mat(ECGdata(locFirstbeat:end, 2)) ./ sampleInt;
    ECGbeatintervals = diff(ECGbeattimes);

end

function [PPGbeattimes, PPGbeatintervals] = readPPG(PPGdata)

    PPGbeattimes = transpose(cell2mat(PPGdata(1, 1)));
    PPGbeatintervals = diff(PPGbeattimes);

end
