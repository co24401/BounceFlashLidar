function [Hpw,Hwp] = compute_homography_matrix(Xw,Xp)
% compute_homography_matrix.m
% Connor Henley
% 9/4/2019
%
% Inputs:
% Xw = Corner points of a plane in world coordinates (2x4)
% Xp = Corner points of plane in pixel coordinates (2x4)
% Order of corner points should be: top-left, top-right, bottom-right,
% bottom-left
%
% Outputs:
% Hpw = Homography matrix, pixel -> world coordiantes (3x3)
% Hwp = Homography matrix, world -> pixel coordiantes (3x3)

% From Tomo's "calibrate.m"

% Construct a transformation equation
A = zeros(8,9);
for ii = 1:4
    A(2*ii-1,:) = [-Xp(1, ii),-Xp(2, ii),-1,0,0,0,Xp(1, ii)*Xw(1, ii),Xp(2, ii)*Xw(1, ii),Xw(1, ii)];
    A(2*ii,:) = [0,0,0,-Xp(1, ii),-Xp(2, ii),-1,Xp(1, ii)*Xw(2, ii),Xp(2, ii)*Xw(2, ii),Xw(2, ii)];
end

% Find the eigen vector with its eigen value = 0
[~,~,V] = svd(A);
Hpw = reshape(V(:,end),[3,3])'; %take the last one
Hwp = inv(Hpw);
end

