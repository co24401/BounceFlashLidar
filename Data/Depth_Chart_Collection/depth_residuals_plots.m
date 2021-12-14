% Boundaries of sets of pixels to use to estimate backplane parameters
umins = [0.1430   -0.0080   -0.0900   -0.1000   -0.0850    0.1080    0.1330];
umaxs = [0.1540    0.1390   -0.0120   -0.0850   -0.0100    0.1440    0.1440];
vmins = [0.0400    0.0800         0   -0.1440   -0.1540   -0.1540   -0.1430];
vmaxs = [0.0887    0.0900    0.0900         0   -0.1430   -0.1430   -0.1270];

inds = false(size(uMap, 1), length(umins));
for ii = 1:size(inds,2)
inds(:, ii) = (uMap > umins(ii)) & (uMap < umaxs(ii)) & (vMap > vmins(ii)) & (vMap < vmaxs(ii));
end
ix_pfit = any(inds, 2);

x_pfit = x_p(ix_pfit,:); % Points on back plane to use for fit
x_pfit = [x_pfit ones(size(x_pfit,1),1)];
[Uf, Sf, Vf] = svd(x_pfit);

backplanef = Vf(:,4);
backplanef = backplanef/norm(backplanef(1:3)); % Normalize first three elements, correspond to plane surface normal

% Compute distances of all points in point cloud from backplane
pdispsf = [x_p ones(size(x_p, 1),1)]*backplanef;

% Plot these distances
figure; imagesc(-1*1000*backplanef(4)*uu, -1*1000*backplanef(4)*vv, reshape(-1*1000*pdispsf, num_u, num_v)'); colorbar;
set(gca, 'XDir', 'reverse', 'YDir', 'normal')
colormap('jet')
caxis(1000*[-.01 .06])

% Bounds of 9 depth squares in uv space
umin_sq = [-0.0700, -0.0740, -0.0740, 0.0070, 0.0070, ...
    0.0020, 0.0820, 0.0820, 0.0750];
umax_sq = [-0.0270, -0.0270, -0.0300, 0.0500, 0.0500, ...
    0.0430, 0.1240, 0.1240, 0.1200];
vmin_sq = [0.0200, -0.0530, -0.1280, 0.0230, -0.0530, ...
    -0.1300, 0.0230, -0.0530, -0.1340];
vmax_sq = [0.0700, -0.0120, -0.0870, 0.0640, -0.0120, ...
    -0.0870, 0.0640, -0.0120, -0.0900];

% Indices of these 9 depth squares
ix_sq = false(size(uMap, 1), length(umin_sq));
for ii = 1:size(ix_sq,2)
ix_sq(:, ii) = (uMap > umin_sq(ii)) & (uMap < umax_sq(ii)) & (vMap > vmin_sq(ii)) & (vMap < vmax_sq(ii));
end

% Array of points corresponding to each depth square
sq_pts = cell(1, size(ix_sq, 2));
for ii = 1:length(sq_pts)
sq_pts{ii} = x_p(ix_sq(:,ii),:);
end

resids = cell(size(sq_pts));
sq_means = zeros(size(sq_pts));
sq_devs = zeros(size(sq_pts));
for ii = 1:length(sq_pts)
resids{ii} = [sq_pts{ii} ones(size(sq_pts{ii},1),1)]*backplanef;
sq_means(ii) = mean(resids{ii});
sq_devs(ii) = std(resids{ii});
end

%%
true_depths = [0 6 12 17 23 29 33 40.5 48];
figure; hold on;
for ii = 1:length(true_depths)
    num_sqpts = length(resids{ii});
    plot(true_depths(ii)*ones(1,num_sqpts), -1000*resids{ii}, '.r')
    plot(true_depths(ii), -1000*sq_means(ii), '_k', 'MarkerSize', 12)
end

all_resids_vs_truth = [];
for ii = 1:length(resids)
    resids_vs_truth = -1000*resids{ii} - true_depths(ii);
    all_resids_vs_truth = [all_resids_vs_truth; resids_vs_truth];
end

stderror_v_truth = sqrt( sum(all_resids_vs_truth.^2) ...
    / (length(all_resids_vs_truth)-1) ); % in mm

%% Individual scan residuals
filedir = 'C:\Users\dolor\Remote-Flash\Experiments\Depth_Chart_Collection_052721\';
filename = 'single_pt_residual_image';

num_illum_pts = 12;
num_squares = length(umin_sq);

singlept_resids = zeros(num_illum_pts, size(r2, 2));

for ii = 1:num_illum_pts

    x_2 = [r2(ii,:)'.*cos(thetaMap) ...
        r2(ii,:)'.*sin(phiMap).*sin(thetaMap) ...
        r2(ii,:)'.*cos(phiMap).*sin(thetaMap)];
    
    singlept_resids(ii,:) = [x_2 ones(size(x_2,1),1)]*backplanef;
    
    % Plot these distances
    figure; imagesc(-1*1000*backplanef(4)*uu, -1*1000*backplanef(4)*vv, reshape(-1*1000*singlept_resids(ii,:), num_u, num_v)'); colorbar;
    set(gca, 'XDir', 'reverse', 'YDir', 'normal')
    colormap('jet')
    caxis(1000*[-.01 .06])
    xlim([-321 431])
    ylim([-431 240])
    xlabel('Horizontal Position (mm)')
    ylabel('Vertical Position (mm)')
    title(['Residuals vs. Back Plane, Point ii = ' num2str(ii)])
    saveas(gcf, [filedir filename '_' num2str(ii)])

end

sq_means_1pt = zeros(num_illum_pts, num_squares);
sq_devs_1pt = zeros(num_illum_pts, num_squares);
stderror_v_truth_1pt = zeros(num_illum_pts,1);
resids_vs_truth_1pt = zeros(num_squares,1);

for ii = 1:num_illum_pts
    resids_vs_truth_1pt = [];
    for jj = 1:num_squares
        resids_1pt = singlept_resids(ii, ix_sq(:, jj));
        sq_means_1pt(ii, jj) = mean(resids_1pt);
        sq_devs_1pt(ii, jj) = std(resids_1pt);
        tmp1 = -1000*resids_1pt - true_depths(jj);
        resids_vs_truth_1pt = [resids_vs_truth_1pt tmp1];
    end
    stderror_v_truth_1pt(ii) = sqrt( sum(resids_vs_truth_1pt.^2) ...
    / (length(resids_vs_truth_1pt)-1) );
end

%% 

rand_ix1 = [6 3 11 7 8 5 1 2 4 9 10 12];  % Random permutation of 12 inds 
% (generated with randperm)

r2_c1 = zeros(num_illum_pts, size(r2, 2));

for ii = 1:num_illum_pts
    cixs = rand_ix1(1:ii);
    r2_c1(ii,:) = sum(R(cixs, :).*S(cixs, :).*r2(cixs, :)./ var_r2(cixs, :), 1, 'omitnan') ./ ...
        sum(R(cixs, :).*S(cixs, :)./var_r2(cixs, :), 1, 'omitnan');
    
end
    
resids_c1 = zeros(num_illum_pts, size(r2_c1, 2));

for ii = 1:num_illum_pts

    x_2 = [r2_c1(ii,:)'.*cos(thetaMap) ...
        r2_c1(ii,:)'.*sin(phiMap).*sin(thetaMap) ...
        r2_c1(ii,:)'.*cos(phiMap).*sin(thetaMap)];
    
    resids_c1(ii,:) = [x_2 ones(size(x_2,1),1)]*backplanef;
    
    % Plot these distances
    figure; imagesc(-1*1000*backplanef(4)*uu, -1*1000*backplanef(4)*vv, reshape(-1*1000*resids_c1(ii,:), num_u, num_v)'); colorbar;
    set(gca, 'XDir', 'reverse', 'YDir', 'normal')
    colormap('jet')
    caxis(1000*[-.01 .06])
    xlim([-321 431])
    ylim([-431 240])
    xlabel('Horizontal Position (mm)')
    ylabel('Vertical Position (mm)')
    title(['Residuals vs. Back Plane, Combined ' num2str(ii) ' pts'])
    saveas(gcf, [filedir 'c1_' filename '_' num2str(ii)])

end

sq_means_c1 = zeros(num_illum_pts, num_squares);
sq_devs_c1 = zeros(num_illum_pts, num_squares);
stderror_v_truth_c1 = zeros(num_illum_pts,1);
stderror_v_mean_c1 = zeros(num_illum_pts,1);

for ii = 1:num_illum_pts
    resids_vs_truth = [];
    resids_vs_mean = [];
    for jj = 1:num_squares
        resids_jj = resids_c1(ii, ix_sq(:, jj));
        sq_means_c1(ii, jj) = mean(resids_jj);
        sq_devs_c1(ii, jj) = std(resids_jj);
        tmp1 = -1000*resids_jj - true_depths(jj);
        resids_vs_truth = [resids_vs_truth tmp1];
        tmp2 = resids_jj-sq_means_c1(ii, jj);
        resids_vs_mean = [resids_vs_mean tmp2];
    end
    stderror_v_truth_c1(ii) = sqrt( sum(resids_vs_truth.^2) ...
    / (length(resids_vs_truth)-1) );
    stderror_v_mean_c1(ii) = sqrt( sum(resids_vs_mean.^2) ...
    / (length(resids_vs_mean)-1) );
end

%%

rand_ix2 = [10 4 5 2 9 7 3 11 6 12 1 8];  % Random permutation of 12 inds 
% (generated with randperm)

r2_c2 = zeros(num_illum_pts, size(r2, 2));

for ii = 1:num_illum_pts
    cixs = rand_ix2(1:ii);
    r2_c2(ii,:) = sum(R(cixs, :).*S(cixs, :).*r2(cixs, :)./ var_r2(cixs, :), 1, 'omitnan') ./ ...
        sum(R(cixs, :).*S(cixs, :)./var_r2(cixs, :), 1, 'omitnan');
    
end
    
resids_c2 = zeros(num_illum_pts, size(r2_c2, 2));

for ii = 1:num_illum_pts

    x_2 = [r2_c2(ii,:)'.*cos(thetaMap) ...
        r2_c2(ii,:)'.*sin(phiMap).*sin(thetaMap) ...
        r2_c2(ii,:)'.*cos(phiMap).*sin(thetaMap)];
    
    resids_c2(ii,:) = [x_2 ones(size(x_2,1),1)]*backplanef;
    
    % Plot these distances
    figure; imagesc(-1*1000*backplanef(4)*uu, -1*1000*backplanef(4)*vv, reshape(-1*1000*resids_c2(ii,:), num_u, num_v)'); colorbar;
    set(gca, 'XDir', 'reverse', 'YDir', 'normal')
    colormap('jet')
    caxis(1000*[-.01 .06])
    xlim([-321 431])
    ylim([-431 240])
    xlabel('Horizontal Position (mm)')
    ylabel('Vertical Position (mm)')
    title(['Residuals vs. Back Plane, Combined ' num2str(ii) ' pts'])
    saveas(gcf, [filedir 'c2_' filename '_' num2str(ii)])

end

sq_means_c2 = zeros(num_illum_pts, num_squares);
sq_devs_c2 = zeros(num_illum_pts, num_squares);
stderror_v_truth_c2 = zeros(num_illum_pts,1);
stderror_v_mean_c2 = zeros(num_illum_pts,1);

for ii = 1:num_illum_pts
    resids_vs_truth = [];
    resids_vs_mean = [];
    for jj = 1:num_squares
        resids_jj = resids_c2(ii, ix_sq(:, jj));
        sq_means_c2(ii, jj) = mean(resids_jj);
        sq_devs_c2(ii, jj) = std(resids_jj);
        tmp1 = -1000*resids_jj - true_depths(jj);
        resids_vs_truth = [resids_vs_truth tmp1];
        tmp2 = resids_jj-sq_means_c2(ii, jj);
        resids_vs_mean = [resids_vs_mean tmp2];
    end
    stderror_v_truth_c2(ii) = sqrt( sum(resids_vs_truth.^2) ...
    / (length(resids_vs_truth)-1) );
    stderror_v_mean_c2(ii) = sqrt( sum(resids_vs_mean.^2) ...
    / (length(resids_vs_mean)-1) );
end