function VXY = galvoCoarsePointAtXYZ(X, Vxshft, Vyshft)
% galvoCoarsePointAtXYZ.m
% Connor Henley
% 11/4/2020
%
% Coarse non-homography based galvo alignment method. Generates command
% voltages that will point at position X in galvo centered coordinates.
% Corrects warping caused by dual-mirror steering.
% Results are usually okay but error can be large (~1cm at 1m depth) at
% large scan angles.  Use a homography to correct further.
%
% Inputs:
% X = Points the you want to aim at with galvo (Nx3)
% Vxshft = Voltage added to compensate for X mirror misalignment (scalar)
% Vyshft = Voltage added to compensate for Y mirror misalignment (scalar)
%
% Outputs:
% VXY = Command voltages that will 

VXY(:, 1) = 0.25*atand(X(:, 1)./sqrt(X(:, 2).^2 + X(:, 3).^2)) + Vxshft;
VXY(:, 2) = 0.25*atand(X(:, 2)./X(:, 3)) + Vyshft;
end

