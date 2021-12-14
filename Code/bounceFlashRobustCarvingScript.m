% bounceFlashRobustCarvingScript.m
% Connor Henley
% 4/29/2021
%
% Adapted from code originally developed for use in:
%
%   C. Henley, T. Maeda, T. Swedish, and R. Raskar, “Imaging behind
%   occluders using two-bounce light,” in Computer Vision – ECCV 2020.
%   Cham: Springer International Publishing, 2020, pp. 573–588.
%
% Please see the above reference for detailed description of algorithm.
%
% Inputs are generated using "twoBounceDetectionProbabilities.m".  
%
% Output "probability_volume" is a voxel grid of occupancy probability
% values.  Output "reconstruction_voxel.volume" is a binary occupancy map,
% and is equivalent to probability_volume > T.

lightShadowFile = 'C:\Users\dolor\Remote-Flash\Experiments\IBO_Collection_051821\two_bounce_detections_Tp9.mat';
load(lightShadowFile, 'illum_pts', 'lit_pts', 'shadow_pts', 'x_p')

%%

%define voxels
outside_voxel.x_lims = [-1 1]; % [-.1 .4]; % [-0.5, 0]; %
outside_voxel.y_lims = [-.6 .6]; %[-.6 -.1];  % [-0.5, 0]; %
outside_voxel.z_lims = [0.5 2.5]; % [1.6 2.1]; % %[1.2 1.7]; 
outside_voxel.num_x = 200; %50; 
outside_voxel.num_y = 120; %50;
outside_voxel.num_z = 200; %50; 
outside_voxel.volume = zeros([outside_voxel.num_x, outside_voxel.num_y, outside_voxel.num_z]);

inside_voxel = outside_voxel;

% reconstruction via carving
disp('Commence robust carving.');
for ii = 1:length(lit_pts)
    [outside_voxel, inside_voxel] = ...
        robust_carving_frame(outside_voxel, inside_voxel, ...
            illum_pts(ii, :), lit_pts{ii}, shadow_pts{ii});
    disp(ii)
end

eta = .05; % Probability occupied voxel is traced to illuminated region (miss probability)
xi = .5; % Probability that an empty voxel is traced to shadow (probability false alarm)
p0 = 0.7; % Prior probability that any voxel is empty
p1 = 0.3; % Prior probability that any voxel is occupied
T = 0.5; % Probability threshold to decide that voxel is occupied

mm = 0:length(lit_pts); 
nn = 0:length(lit_pts);
[M, N] = meshgrid(mm, nn);
PO = p1*(eta.^M).*((1-eta).^N)./(p0*((1-xi).^M).*(xi.^N) + p1*(eta.^M).*((1-eta).^N));

% mm = 0:(2*length(lit_pts)); 
% nn = 0:(2*length(lit_pts));
% [M, N] = meshgrid(mm, nn);
% PO = p1*(eta.^M).*((1-eta).^N)./(p0*((1-xi).^M).*(xi.^N) + p1*(eta.^M).*((1-eta).^N));

figure;
[reconstruction_voxel.volume, probability_volume] = ...
    visualize_probablistic(inside_voxel, outside_voxel ,PO, T, true);
