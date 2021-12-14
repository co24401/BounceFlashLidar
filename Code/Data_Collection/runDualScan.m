% runDualScan.m
% Connor Henley
% 2/26/2021
%
% Script to run data collection and output .csv files that contain raw
% photon count data from Hydraharp.

%% Galvo and Hydraharp setup

NDF = true; % Set to true for attenuated data collection, false for 
% regular data collection

% Set up galvos

cscanfile = 'C:\Users\dolor\Remote-Flash\Experiments\TCI_large_scene_16pts\detector_scan_200x200.mat';
load(cscanfile, 'V_UV')
V_UVc = V_UV;
clear V_UV

lscanfile = 'C:\Users\dolor\Remote-Flash\ExperimentsTCI_large_scene_16pts\laser_scan_16pts.mat';
load(lscanfile, 'V_UV')
V_UVl = V_UV;
clear V_UV

sc = daq.createSession('ni');
addAnalogOutputChannel(sc, 'Dev1', 'ao0', 'Voltage')
addAnalogOutputChannel(sc, 'Dev1', 'ao1', 'Voltage')
addDigitalChannel(sc, 'Dev1', 'port0/line2', 'OutputOnly')
outputSingleScan(sc, [V_UVc(1, 1), V_UVc(1, 2), 0])

sl = daq.createSession('ni');
addAnalogOutputChannel(sl, 'Dev2', 'ao0', 'Voltage')
addAnalogOutputChannel(sl, 'Dev2', 'ao1', 'Voltage')
outputSingleScan(sl, [V_UVl(1, 1), V_UVl(1, 2)])

dwellTime = 0.1;
Tacq = dwellTime*size(V_UVc, 1)+20;
disp(['Hydraharp measurement interval: ' num2str(Tacq)]);
disp(['Combined time for all scans: ' num2str(size(V_UVl, 1)*Tacq)]);
sc.Rate = 2/dwellTime;

% Set of NI commands to execute detector galvo scan sequenc with NI card
cmds = [kron(V_UVc, [1; 1]), repmat([1; 0], size(V_UVc, 1), 1)];
cmds = [cmds; [V_UVc(end, 1), V_UVc(end, 2), 1] ];

% Hydraharp acquisition settings
Binning       = 3;  % Resolution is (base resolution)^(Binning) (eg 2^3 = 8 ps)
SyncDiv       = 4;  % See HH manual for description of remaining parameters
SyncCFDZeroX  = 10;
SyncCFDLevel  = 400;
InputCFDZeroX = 10;
InputCFDLevel = 300;

bin_size = 8E-12;

datapath = 'C:\Users\dolor\Remote-Flash\Experiments\mirror_scan_100721\';

%% Scan and acquire

for ii = 2:size(V_UVl, 1)
    
    disp(['Scan to illumination point no. ' num2str(ii) '.'])
    
    if NDF
        filename = ['spot_' num2str(ii) '_NDF.out'];
    else
        filename = ['spot_' num2str(ii) '.out'];
    end
    
    % Next laser spot
    outputSingleScan(sl, V_UVl(ii, :))
    % Reset detector galvo scan to starting position
    outputSingleScan(sc, [V_UVc(1, 1), V_UVc(1, 2), 0])
    % Queue detector galvo voltage commands
    queueOutputData(sc, cmds)
    pause(1)
    
    disp('Begin data collection.')
    
    syncrate = HHscanAndRecordT3([datapath filename], sc, ...
        1000*Tacq, Binning, ...
        SyncDiv, SyncCFDZeroX, SyncCFDLevel, ...
        InputCFDZeroX, InputCFDLevel);
    
    pause(5)
    
    disp(['Converting file ' num2str(ii) ' to .csv.'])
    
    if NDF
        filename = ['spot_' num2str(ii) '_NDF.out'];
    else
        filename = ['spot_' num2str(ii) '.out'];
    end
    
    global_resolution = 1/syncrate;
    
    Read_T3OUT_csv([datapath filename], bin_size, global_resolution)
    
end

disp('Data acquisition complete')