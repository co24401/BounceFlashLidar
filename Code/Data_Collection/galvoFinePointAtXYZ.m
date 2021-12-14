function VXY = galvoFinePointAtXYZ(Xa, Hau, Vxshft, Vyshft)
% galvoCoarsePointAtXYZ.m
% Connor Henley
% 11/4/2020
%
% Generates command voltages that will point at position Xa in galvo 
% centered coordinates.  Does this by first using a homography (generated 
% previously in calibration step by user) to convert from true (aligned)
% coordinates to unaligned galvo-centered coordinates (if you tried to
% point galvos assuming no misalignment, they would point at Xu, not Xa).
% These unaligned coordinates can be converted to command voltages with an
% analytical formula that relates mirror tilt angles to pointing direction.
%
% Inputs:
% Xa = Points the you want to aim at with galvo (Nx3)
% Hau = Homography that transforms from aligned to unaligned frame (3x3)
% Vxshft = Voltage added to compensate for X mirror misalignment (scalar)
% Vyshft = Voltage added to compensate for Y mirror misalignment (scalar)
%
% Outputs:
% VXY = Command voltages that will result in galvo pointing at Xa (Nx2)

homog = @(x, H) H*x./(H(3, :)*x);

Xu = homog(Xa', Hau)';

VXY = galvoCoarsePointAtXYZ(Xu, Vxshft, Vyshft);

end

