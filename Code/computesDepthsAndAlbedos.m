% computeDepthsAndAlbedos.m
% Connor Henley
% 4/16/2021
%
% Using inputs from "batchProcessRemoteFlash.m" (typically saved as 
% "params.mat"), generates range ("r2_est") and albedo (rho2_est) maps of 
% scene, and also combines range map with galvo scan angle information to 
% save a point cloud representation of scene (x_p).  Saves outputs in a
% .mat file (typically named "results.mat").
%
% Outputs from this script are used as inputs for 
% "twoBounceDetectionProbabilities.m".

filedir = 'C:\Users\dolor\Remote-Flash\Experiments\Depth_Chart_Collection_052721\';
filename = 'results.mat'; % Filename of saved output .mat
paramfilename = 'params.mat';
scanfile = 'C:\Users\dolor\Remote-Flash\Experiments\Depth_Chart_Collection_052721\detector_scan_135x135.mat';

load([filedir paramfilename])
load(scanfile, 'num_u', 'num_v', 'u_lims', 'v_lims', 'thetaMap', 'phiMap')

t0 = 14.2428e-9; % Time-of-arrival associated with zero time-of-flight
bin_width = 8E-12; % picoseconds

c = 299792458; % speed of light (m/s)
s = .257; % baseline separation between galvos
OD = 100; % NDF filter attenuation factor
dwell_time = 0.05;  % per-pixel dwell time (sec)

B1 = dwell_time*1.514e9/9;%dwell_time*1.514e9;
B2 = dwell_time*4.96e4;

halo_bins = 10; % If return within halo_bins of 1B return, assume is lens flare and remove
w2_min = 2; % Min and max acceptable (pulse-width^2) in timing bin units
w2_max = 185; % Used to remove 3+Bounce returns misclassified as 2B

ixs1 = zeros(size(Evals,1), 1);
E1 = zeros(size(Evals,1), 1);
var_E1 = zeros(size(Evals,1), 1);
t1 = zeros(size(Evals,1), 1);
var_t1 = zeros(size(Evals,1), 1);
r1 = zeros(size(Evals,1), 1);
var_r1 = zeros(size(Evals,1), 1);
theta1 = zeros(size(Evals,1), 1);
phi1 = zeros(size(Evals,1), 1);
dr1_dt1 = zeros(size(Evals,1), 1);

t2 = bin_width*pkbins - t0; % Fitted peak bin converted to time-of-flight (in sec)
var_t2 = bin_width*bin_width*muvars; % Variance of t2 estimates
E2 = zeros(size(Evals));
var_E2 = zeros(size(Evals));
r2 = zeros(size(Evals));
var_r2 = zeros(size(Evals));

%% Process single bounce returns

disp('Processing single-bounce returns.')
for ii = 1:size(Evals_NDF, 1)
    
    % Pixel with max Eval assumed to be laser spot
    [mx, ix] = max(Evals_NDF(ii, :));
    [row, col] = ind2sub([num_u num_v], ix);
    ixs1(ii) = ix;
    
    theta1(ii) = thetaMap(ix);
    phi1(ii) = phiMap(ix);
    
    % Assume beam energy in a 5x5 window surrounding max. Eval estimate.
    rows = row + repmat([-2 -1 0 1 2]', 1, 5);
    cols = col + repmat([-2 -1 0 1 2], 5, 1);
    inbounds = rows > 0 & rows <= num_u & cols > 0 & cols <= num_v;
    rows = rows(inbounds);
    cols = cols(inbounds);
    
    beam_ixs = sub2ind([num_u num_v], rows, cols);
            
    E1(ii) = sum(Evals_NDF(ii, beam_ixs), 'omitnan');
    var_E1(ii) = sum(Evars_NDF(ii, beam_ixs), 'omitnan');
    pk1 = sum(S_NDF(ii, beam_ixs).*pkbins_NDF(ii, beam_ixs)./ ...
        muvars_NDF(ii, beam_ixs), 'omitnan') ./ ...
        sum(S_NDF(ii, beam_ixs) ./ muvars_NDF(ii, beam_ixs), 'omitnan');
    t1(ii) = bin_width*pk1 - t0;
    var_pk1 = ...
        1 / sum(S_NDF(ii, beam_ixs) ./ muvars_NDF(ii, beam_ixs), 'omitnan');
    var_t1(ii) = bin_width*bin_width*var_pk1;
    
    cos_theta1c = cos(theta1(ii));
    denom = c*t1(ii) - s*cos_theta1c;
    
    r1(ii) = 0.5*(c^2 * t1(ii)^2 - s^2)/denom;
    dr1_dt1(ii) = (c/2)*( 1 + (s^2)*(1 - cos_theta1c^2)/(denom^2) );
    var_r1(ii) = (dr1_dt1(ii)^2)*(var_t1(ii)^2);
    
end

no_1B_halo = pkbins > (halo_bins + (t1+t0)/bin_width);  % Filter that removes lens flare artifacts
no_3B = w2s > w2_min & w2s < w2_max; % Filter that supresses likely 3B estimates
R = no_1B_halo & no_3B;

%% Process two-bounce returns

disp('Processing two-bounce returns.')
for ii = 1:size(Evals, 1)
   
    dt12 = t2(ii, :) - t1(ii);
    cos_delta = cos(theta1(ii))*cos(thetaMap) + ...
        sin(theta1(ii))*sin(thetaMap).*cos(phi1(ii)-phiMap);
    
    a = (1 - cos_delta'.^2)*(r1(ii)^2)/(c^2);
    b = dt12 + (1-cos_delta')*r1(ii)/c;
    
    r2(ii, :) = (c/2)*( dt12.*(dt12 + 2*r1(ii)/c) ) ./ b;
    dr2_dt1 = (c/2) * ( ((1 + cos_delta').*(dt12.^2).* dr1_dt1(ii)/c - a) ...
        ./ (b.^2) - 1);
    dr2_dt2 = (c/2) * (1 + a ./ (b.^2));
    
    var_r2(ii, :) = (dr2_dt1.^2).*var_t1(ii) + (dr2_dt2.^2).*var_t2(ii, :);
    
end

r2_est = sum(R.*S.*r2./ var_r2, 1, 'omitnan') ./ ...
    sum(R.*S./var_r2, 1, 'omitnan');  % Is this going to miscount when e.g. numerator is nan and denominator isn't?
var_r2_est = 1 ./ sum(R.*S./var_r2, 1, 'omitnan');

%% Point cloud and normals

% Create point cloud
x_p = [r2_est'.*cos(thetaMap) ...
    r2_est'.*sin(phiMap).*sin(thetaMap) ...
    r2_est'.*cos(phiMap).*sin(thetaMap)];

figure; scatter3(x_p(:, 1), x_p(:, 2), x_p(:, 3), 6, r2_est')
xlim([-1 1])
ylim([-1 1])
zlim([.5 2.5])
xlabel('X')
ylabel('Y')
zlabel('Z')
caxis([1 3])

pc = pointCloud(x_p);
 
n_p = pcnormals(pc, 20);

x_c = [0,0,0];
for ii = 1:size(x_p, 1)
   p1 = x_c - x_p(ii, :);
   p2 = n_p(ii, :);
   % Flip the normal vector if it is not pointing towards the sensor.
   angle = atan2(norm(cross(p1,p2)),p1*p2');
   if angle > pi/2 || angle < -pi/2
       n_p(ii, :) = -n_p(ii, :);
   end
end

%% Albedo estimation

rho1 = zeros(size(Evals,1), 1);
var_rho1 = zeros(size(Evals,1), 1);
n1 = zeros(size(Evals,1), 3);

rho2 = zeros(size(Evals));
var_rho2 = zeros(size(Evals));

n_c = [0 0 1];
min_cos = .01;
max_hn = 1;

for ii = 1:size(Evals_NDF, 1)
    
    [mx, ix] = max(Evals_NDF(ii, :));
   
    n1 = n_p(ix, :);
    
    r1c = norm(x_c-x_p(ix, :));
    w1c = (x_c-x_p(ix, :))/r1c;
    %a1 = OD*(w1c*n1')*(-w1c*n_c')/(pi*r1c^2);
    a1 = (1/OD)*B1*(w1c*n1')*(-w1c*n_c')/(pi*r1c^2);
    
    rho1(ii) = E1(ii)/a1;
    var_rho1(ii) = var_E1(ii)/a1^2;
    
    r12 = vecnorm(x_p-x_p(ix, :), 2, 2);
    w12 = (x_p-x_p(ix, :))./r12;
    
    
    cos12 = sum(w12.*n_p(ix, :), 2);
    cos12(cos12<min_cos) = NaN;  % If line-of-sight between 1 and 2 passes behind 
                      % surface, we know it's occluded.  Set term to zero.
    cos21 = dot(-w12, n_p, 2);
    cos21(cos21<min_cos) = NaN; 
    r2c = vecnorm(x_p-x_c, 2, 2);
    wc2 = (x_p - x_c)./r2c;
    cosc2 = sum(wc2.*n_c, 2);
    cosc2(cosc2<min_cos) = NaN; 
    
    h2 = -(wc2 + w12)./vecnorm(wc2 + w12, 2, 2);
    h2n2 = sum(h2.*n_p, 2);
    specular = ones(size(cos21));
    specular(h2n2>max_hn) = NaN;
    
%    a12 = (dPix/pi^2)*cos12.*cos21.*cosc2./r12.^2;
    a12 = B2*(1/pi^2)*cos12.*cos21.*cosc2.*specular./r12.^2;
    
    rho2(ii, :) = Evals(ii, :)./(rho1(ii)*a12');
    var_rho2(ii, :) = Evars(ii, :)./(rho1(ii)*a12').^2 + ...
                var_rho1(ii)*(Evals(ii, :).^2)./ ((rho1(ii)^4)*(a12'.^2));
    
end

rho2_est = sum(R.*S.*rho2./var_rho2, 1, 'omitnan') ./ sum(R.*S.*sign(rho2)./var_rho2, 1, 'omitnan');
var_rho2_est = 1 ./ sum(R.*S./var_rho2, 1, 'omitnan');
rho2_avg = sum(R.*S.*rho2, 1, 'omitnan') ./ sum(R.*S.*sign(rho2), 1, 'omitnan');

figure; imagesc(reshape(rho2_est, num_u, num_v)'); colorbar;
set(gca, 'XDir', 'normal', 'YDir', 'normal')
caxis([0 1.5])
title(['max h2n2 = ' num2str(max_hn)])

figure; imagesc(reshape(rho2_avg, num_u, num_v)'); colorbar;
set(gca, 'XDir', 'normal', 'YDir', 'normal')
caxis([0 1.5])
title(['naive average, max h2n2 = ' num2str(max_hn)])

%% Save results

save([filedir filename], 'r2_est', 'rho2_est', 'r1', 'r2', 'var_r1', 'var_r2', ...
'var_r2_est', 'rho1', 'rho2', 'var_rho1', 'var_rho2', 'var_rho2_est', 'n_c', ...
'n_p', 'R', 'ixs1', 'x_c')