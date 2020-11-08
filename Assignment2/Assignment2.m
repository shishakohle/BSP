%% Biosignal Processing - Assingment 2
% Group H: Alexander Neumayer, Ingo Weigel
% 2020-11-08
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

[beatTimesAmplitudesPQIs, pulseWaveTemplate] = PQI(rawPPGsignal, ...
    samplingRate);
