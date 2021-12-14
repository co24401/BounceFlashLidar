function [Hpw, Hwp] = galvoHomography(Vg, Xw, Vxshft, Vyshft)
% galvoHomography.m
% Connor Henley
% 10/3/2019
%
% Inputs:
% Vg = Voltage commands required to aim laser at four corners of plane (2x4)
% Xw = Corner points of a plane in world coordinates (2x4)
% Vyshft = Voltage added to compensate for mirror misalignment
%
% Outputs:
% Hpw = Homography matrix, pixel -> world coordiantes (3x3)
% Hwp = Homography matrix, world -> pixel coordiantes (3x3)

    Xp = [tand(4*(Vg(1, :)-Vxshft))./cosd(4*(Vg(2, :)-Vyshft)); tand(4*(Vg(2, :)-Vyshft))];

    [Hpw, Hwp] = compute_homography_matrix(Xw, Xp);

end

