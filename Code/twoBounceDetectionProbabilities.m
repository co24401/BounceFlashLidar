% twoBounceDetectionProbabilities.m
% Connor Henley
% 4/29/2021
%
% Using inputs from "computeDepthsAndAlbedos.m" (typically saved as 
% "results.mat"), generates gated two-bounce images (Q2B_vals), map of
% expected two-bounce energies based on scattering point position and
% albedo, and assuming no shadow (E_Q2B_vals).  Map of probabilities that
% each pixel contains a two-bounce signal (i.e. it isn't in shadow, these 
% probabilities are saved in Px1y). These structures all have dimensions of
% (num. illumination spots) x (num. of pixels).
%
% For each illumination point position, also generates point cloud of
% not-shadowed points (lit_pts), and shadowed points (shadow_pts), and
% saves laser spot positions (illum_pts).
%
% Outputs (typically saved as "two_bounce_detections.mat") are used by 
% "bounceFlashRobustCarvingScript.m".

% File directory with the data cube
filedir = 'C:\Users\dolor\Remote-Flash\Experiments\IBO_Collection_051821\';
lightShadowFile = 'two_bounce_detections_Tp9.mat';

% File directory with range and albedo maps
mapFile = 'C:\Users\dolor\Remote-Flash\Experiments\IBO_Collection_051821\results.mat';
load(mapFile, 'r2_est', 'var_r2_est', 'rho2_est', 'var_rho2_est', 'r1', 'var_r1', 'rho1', 'var_rho1', 'ixs1', 'n_p');

% File directory with galvo scan parameters
scanFile = 'C:\Users\dolor\Remote-Flash\Experiments\IBO_Collection_051821\detector_scan_200x200.mat';
load(scanFile, 'thetaMap', 'phiMap')

NDF = false;

t0 = 14.2428e-9; % Time-of-arrival associated with zero time-of-flight
bin_width = 8E-12; % picoseconds

bin_radius = 30; 
noise_gate_start = 1; 
noise_gate_stop = 1000;

c = 299792458; % speed of light (m/s)
s = .257; % baseline separation between galvos
OD = 100; % NDF filter attenuation factor
dwell_time = 0.1;  % per-pixel dwell time (sec)

B1 = dwell_time*1.514e9/9;%dwell_time*1.514e9;
B2 = dwell_time*4.96e4;

if NDF
    B2 = B2/OD;
end

min_cos = 0.01;
x_c = [0 0 0];
n_c = [0 0 1];

% Generate point cloud from depth map
x_p = [r2_est'.*cos(thetaMap) ...
    r2_est'.*sin(phiMap).*sin(thetaMap) ...
    r2_est'.*cos(phiMap).*sin(thetaMap)];

illum_pts = zeros(length(r1), 3);
lit_pts = cell(length(r1), 1);
shadow_pts = cell(length(r1), 1);
Q2B_vals = zeros(length(r1), length(r2_est));
E_Q2B_vals = zeros(length(r1), length(r2_est));

% Shadow thresholding parameters
  
p0 = 0.5; % Prior probability that 2B signal is not present
p1 = 0.5; % Prior probability that 2B signal is present

T = 0.9; %0.5; % 2B detection probability threshold

Px1y = zeros(length(r1), length(r2_est));

for ii = 1:length(r1)
    
    if NDF
        filename = ['spot_' num2str(ii) '_NDF'];
    else
        filename = ['spot_' num2str(ii)];
    end
    
    disp(['Laser position' num2str(ii)])
    
    load([filedir filename '.mat'], 'dataCube')

    ix1 = ixs1(ii);
    
    % Calculate expected 2B return time

    theta1c = thetaMap(ix1);
    phi1c = phiMap(ix1);
    r1c = r1(ii);
    
    r1l = sqrt(r1c^2 + s^2 - 2*r1c*s*cos(theta1c));
    
    cos_delta = cos(theta1c)*cos(thetaMap) + ...
        sin(theta1c)*sin(thetaMap).*cos(phi1c-phiMap);
    
    r12s = sqrt(r1c^2 + r2_est.^2 - 2*r1c*r2_est.*cos_delta');
    
    t2B = (r1l + r12s + r2_est) / c;
    
    bin_2B = (t2B + t0) / bin_width;
    
    % Estimate background
    b = mean(dataCube(:, noise_gate_start:noise_gate_stop), 2);
    b = max(b, 1/(noise_gate_stop - noise_gate_start + 1)); % If estimated 
    % background is zero, any counts will trigger P = 1 2B detection.
    % Avoid this by setting b to minimum possible detectable value if
    % observed to be zero
    
    Q2B = nan(size(b));
    
    for jj = 1:size(dataCube, 1)
        cen = round(bin_2B(jj));
        if ~isnan(cen)
            start = max(cen-bin_radius, 1);
            stop = min(cen+bin_radius, size(dataCube, 2));
            Q2B(jj) = sum(dataCube(jj, start:stop));
        end
    end
    
    Q2B_vals(ii, :) = Q2B;
    
    x_1 = r1c * [cos(theta1c) sin(phi1c)*sin(theta1c) cos(phi1c)*sin(theta1c)];
    illum_pts(ii, :) = x_1;
    
    n_1 = n_p(ix1, :);
    
    w12 = (x_p-x_1)./r12s';
    
    cos12 = sum(w12.*n_1, 2);
    cos12(cos12<min_cos) = NaN;  % If line-of-sight between 1 and 2 passes behind
    % surface, we know it's occluded.  Set term to zero.
    cos21 = dot(-w12, n_p, 2);
    cos21(cos21<min_cos) = NaN;
    wc2 = (x_p - x_c)./r2_est';
    cosc2 = sum(wc2.*n_c, 2);
    cosc2(cosc2<min_cos) = NaN;
    
    G12 = (1/pi^2)*cos12.*cos21.*cosc2./r12s'.^2;
    
    if isnan(rho2_est(ix1))
        rho1_est = rho1;
    else
        rho1_est = rho2_est(ix1);
    end
    
    E_Q2B = B2*rho1_est*G12.*rho2_est'; % Expected two-bounce energy
    R = ((r12s/c)/bin_width) < bin_radius ;  % To reject pixels that likely
    % contain significant 1B halo returns
    
    E_Q2B_vals(ii, :) = E_Q2B;
    
    % Probability based thresholding.
    
    lam0 = 2*bin_radius*b;
    lam1 = E_Q2B + 2*bin_radius*b;
    
    q = log(p0) - log(p1) + Q2B.*log(lam0) - Q2B.*log(lam1) - lam0 + lam1;   
    Px1y(ii, :) = 1 ./ (1 + exp(q));
    
    lit_pts{ii} = x_p(Px1y(ii, :) > T, :);
    shadow_pts{ii} = x_p(Px1y(ii, :) < T, :);
end

has_2B = Px1y > T;

% Add on an extra illum. pt at the camera.
illum_pts(end+1, :) = x_c;
lit_pts{end+1} = x_p;
shadow_pts{end+1}=[];

save([filedir lightShadowFile], 'Q2B_vals', 'E_Q2B_vals', 'Px1y', 'has_2B', ...
   'illum_pts', 'lit_pts', 'shadow_pts', 'T', 'p0', 'p1', 'x_p')