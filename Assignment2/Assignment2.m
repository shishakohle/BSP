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
filepath = edfSubfolder + "\" + edfFilename;
% uncomment the next line for Linux distributions / Mac OS X
% filepath = edfSubfolder + "/" + edfFilename;

clear edfSubfolder;
clear edfFilename;

%% aquire signal
[hdr, record] = edfread(filepath); %% TODO record needed?
[hdr, rawPPGsignal] = edfread(filepath,'targetSignals','wrist_ppg');
samplingRate = hdr.samples/hdr.duration;

pqi = PQI(rawPPGsignal, samplingRate, record);
