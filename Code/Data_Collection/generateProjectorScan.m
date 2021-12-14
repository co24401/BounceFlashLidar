function [V_UV, UV] = generateProjectorScan(Hau, num_u, num_v, u_lims, v_lims, Vxshft, Vyshft)
% generateProjectorScan.m
% Connor Henley
% 11/4/2019
%
% Inputs:
% Hau = Homography matrix, aligned -> unaligned frame (3x3)
% num_u = Number of columns in grid scan (scalar int)
% num_v = Number of rows in grid scan (scalar int)
% u_lims = Vector containing min and max desired u values (e.g. [umin umax]) (1x2)
% v_lims = Vector containing min and max desired v values (1x2)
% Vxshft = Voltage added to compensate for X-mirror misalignment
% Vyshft = Voltage added to compensate for Y-mirror misalignment
%
% Outputs:
% V_UV = Array of voltage commands to be sent to galvo.  Scan will start in
%   top left corner and snake down to bottom left or right corner,
%   alternating scan directions at each row
% UV1 = UV pixel coordinates corresponding to each voltage
%   command

uu = linspace(u_lims(1), u_lims(2), num_u);
vv = linspace(v_lims(1), v_lims(2), num_v);

UU = repmat(uu', 1, length(vv));
VV = repmat(vv, length(uu), 1);

for ii = 1:size(UU, 2)
    if mod(ii, 2)
        UU(:, ii) = flipud(UU(:, ii));
        VV(:, ii) = flipud(VV(:, ii));
    end
end

UV = [UU(:) VV(:) ones(size(UU(:)))];

V_UV = galvoFinePointAtXYZ(UV, Hau, Vxshft, Vyshft);

if abs(max(V_UV(:))) > 10
    disp('Warning: Maximum permissible voltage exceeded.')
end

end

