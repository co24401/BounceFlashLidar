% calibrateGalvoAsProjector.m
% Connor Henley
% 10/3/2019
%
% Script used to calibrate galvo mirrors.  Requires determining command 
% voltages that prompt a galvo system to point at four points on a plane
% with known position (in world coordinates). See subfunction comments for
% descriptions of specific inputs and outputs.

%%
folder = 'C:\Users\dolor\Remote-Flash\Experiments\mirror_scan_100621\';
filename = 'laser_scan_10x4.mat';

%% Corner point galvo votlage measurements
% Laser galvo cal points
Vg = [-5.09, -1.2; 5.62, -1.13; 5.25, 7.08; -4.72, 7.02]';
% % Detector galvo cal points
% Vg = [-2.592, -0.591; 7.794, -0.591; 7.592, 7.490; -2.088, 7.490]';

% Laser galvo shifts
Vxshft = 0.8;  
Vyshft = 0.97; % Offset to approximately account for galvo misalignment
% Detector galvo shifts
% Vxshft = 0.53;  
% Vyshft = 1.5; % Offset to approximately account for galvo misalignment

%% Compute homographies that map projector coordiantes to each observation wall

% Generate homography

% laser galvo coordinate corner points (galvo centered, cm)
A = [-56.5, -19, 128.5]; % Bottom right corner
B = [45.2, -19, 128.5]; % Bottom left corner
C = [45.2, 57.0, 128.5]; % Top left corner
D = [-56.5, 57.0, 128.5]; % Top right corner
Xa = [A; B; C; D];

% % Detector galvo coordinate corner points (galvo centered, cm)
% A = [-30.8, -18.9, 128.5]; % Bottom right corner
% B = [70.9, -18.9, 128.5]; % Bottom left corner
% C = [70.9, 57.1, 128.5]; % Top left corner
% D = [-30.8, 57.1, 128.5]; % Top right corner
% Xa = [A; B; C; D];

[Hua, Hau] = galvoHomography(Vg, (Xa(:, 1:2)./Xa(:, 3))', Vxshft, Vyshft);

%% Define scan grids

num_u = 10;
num_v = 4;
u_lims = [-.695 .155];
v_lims = [-.25 0];

[V_UV, UV] = generateProjectorScan(Hau, num_u, num_v, u_lims, v_lims, Vxshft, Vyshft);

%% Calculate phi, theta for future processing
uMap = UV(:, 1);
vMap = UV(:, 2);

% Flip every other row left-right
for ii = 1:num_v
    if rem(ii, 2) == 0
        flip_inds = (1:num_u) + (ii-1)*num_u;
        uMap(flip_inds) = flip(uMap(flip_inds), 1);
    end
end

thetaMap = acot(uMap ./ sqrt(1 + vMap.^2));
thetaMap(thetaMap<0) = thetaMap(thetaMap<0)+pi;
phiMap = atan(vMap);

%%
save([folder filename], 'V_UV', 'UV', 'uMap', 'vMap', 'thetaMap', 'phiMap', 'Hua', 'Hau', 'Vxshft', 'Vyshft', 'Vg', 'Xa', 'num_u', 'num_v', 'u_lims', 'v_lims')