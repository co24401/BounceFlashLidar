% batchProcessRemoteFlash.m
% Connor Henley
% 4/26/2021
%
% Loop through raw photon count data to generate per-pixel statistics such
% as peak bin position, photon counts associated with peak return, pulse
% width, and whether or not a pixel is in shadow.
%
% Outputs from this scripts ("params.mat") used as input to depth and
% albedo estimation script ("computeDepthsAndAlbedos.m")

filedir = 'C:\Users\dolor\Remote-Flash\Experiments\Depth_Chart_Collection_052721\';
scanfile = 'C:\Users\dolor\Remote-Flash\Experiments\Depth_Chart_Collection_052721\detector_scan_135x135.mat';
replicafile = 'C:\Users\dolor\Remote-Flash\Experiments\Replica Collection 042821\replica.mat';

paramfilename = 'params.mat';  % File to be saved as output

% Parameters for removing blur caused by galvo motion
dwell = .05; % seconds
lead_buffer = 0.0005; % seconds
trail_buffer = 0; % seconds
prf = 20E6; % Hz

load(scanfile, 'num_u', 'num_v')
num_pix = num_u*num_v;
num_bins = 8192;
num_pts = 16;

load(replicafile, 'mu_fit_LUT', 'pkbias_LUT')

savefigs = true;
savedata = true;

bin_radius = 30; 
noise_gate_start = 1; 
noise_gate_stop = 1000;
min_pk_cts_threshold = 0.8;%2;

pkbins = zeros(num_pts, num_pix);
mubins = zeros(num_pts, num_pix);
muvars = zeros(num_pts, num_pix);
Evals = zeros(num_pts, num_pix);
Evars = zeros(num_pts, num_pix);
w2s = zeros(num_pts, num_pix);
S = false(num_pts, num_pix);

pkbins_NDF = zeros(num_pts, num_pix);
mubins_NDF = zeros(num_pts, num_pix);
muvars_NDF = zeros(num_pts, num_pix);
Evals_NDF = zeros(num_pts, num_pix);
Evars_NDF = zeros(num_pts, num_pix);
w2s_NDF = zeros(num_pts, num_pix);
S_NDF = false(num_pts, num_pix);

warning('off', 'MATLAB:nearlySingularMatrix')

disp('Processing NDF collects')
%%
for ii = 1:num_pts
    
    filename = ['spot_' num2str(ii) '_NDF'];
    
    disp(['Laser position' num2str(ii)])
    
    [dataCube, ~, ~, ~] = processRawCounts( ...
       filename, filedir, scanfile, savedata, savefigs, ...
       dwell, lead_buffer, trail_buffer, prf, num_bins);
    
    [pkbins_NDF(ii, :), mubins_NDF(ii, :), muvars_NDF(ii, :), Evals_NDF(ii, :), Evars_NDF(ii, :), ...
        w2s_NDF(ii, :), S_NDF(ii, :)] = datacubeStats(dataCube, ...
        noise_gate_start, noise_gate_stop, ...
        min_pk_cts_threshold, mu_fit_LUT, pkbias_LUT);
    
end

disp('Processing no-NDF collects')

for ii = 1:num_pts
    
    filename = ['spot_' num2str(ii)];
    
    disp(['Laser position' num2str(ii)])
    
    [dataCube, ~, ~, ~] = processRawCounts( ...
       filename, filedir, scanfile, savedata, savefigs, ...
       dwell, lead_buffer, trail_buffer, prf, num_bins);
    
    [pkbins(ii, :), mubins(ii, :), muvars(ii, :), Evals(ii, :), Evars(ii, :), ...
        w2s(ii, :), S(ii, :)] = datacubeStats(dataCube, ...
        noise_gate_start, noise_gate_stop, min_pk_cts_threshold, ...
        mu_fit_LUT, pkbias_LUT);
    
end

warning('on', 'MATLAB:nearlySingularMatrix')

save([filedir paramfilename], 'pkbins', 'pkbins_NDF', ...
    'mubins', 'muvars', 'Evals', 'Evars', 'w2s', 'S', ...
    'mubins_NDF','muvars_NDF','Evals_NDF', 'Evars_NDF', 'w2s_NDF', 'S_NDF')
