function [p_pk, p_mu, var_p_mu, E, var_E, p_w2, S] = ...
    datacubeStats(dataCube, noise_gate_start, noise_gate_stop, ...
    threshold_absolute, mu_fit_LUT, pkbias_LUT)
% dataCubeStats.m
% Connor Henley
% 4/26/2021
%
%  Reads in dataCube and computes signal statistics for each histogram.

% Inputs:
% dataCube = Counts binned by pixel and timing position
% noise_gate_start = First bin in timing gate used to estimate background
%   noise level in a histogram
% noise_gate_stop = Last bin in timing gate used to estimate background
%   noise level in a histogram
% threshold_absolute = Peak of histogram convolved with Gaussian matched
%   filter must be above this value, otherwise the pixel that corresponds
%   to that histogrma will be declared "in shadow"
% mu_fit_LUT = Lookup table for correcting estimation bias due to finite 
%   fitting window.  Corrects biases in estimated mean timing bin for 
%   two-bounce returns. 
% pkbias_LUT = = Lookup table for correcting estimation bias due to finite 
%   fitting window.  Corrects biases in estimated peak timing bin for 
%   two-bounce returns. 

% Outputs:
% p_pk = Pixel-wise map of estimated peak bins
% p_mu = Pixel-wise map of estimated mean bin (within fitting window)
% var_p_mu = Estimated uncertainty (variance) in p_mu estimates
% E = Pixel-wise map of estimated two-bounce energy
% var_E = Estimated uncertainty (variance) in two-bounce energy estimates
% p_w2 = Pixel-wise map of estimated two-bounce return pulsewidth
% S = Pixel-wise shadow classification (0 = in shadow, 1 = not in shadow)

num_noise_bins = noise_gate_stop - noise_gate_start + 1;

% Gaussian used for convolution / window finding
w2_f = 51.1;
filter_halfwidth = 50;
filter_halfwidth = ceil(filter_halfwidth);
tt = -filter_halfwidth:filter_halfwidth;
ff = (1 / sqrt(2*pi*w2_f))*exp(-(tt.^2) ./ (2*w2_f));

conv2_datacube = conv2(dataCube, ff, 'same');

[pkvals, pkbins] = max(conv2_datacube, [], 2, 'omitnan');

b = mean(dataCube(:, noise_gate_start:noise_gate_stop), 2);
E = nan(size(b));
p_mu = nan(size(b));
p_pk = nan(size(b));
p_w2 = nan(size(b));

var_p_mu = nan(size(b));
var_E = nan(size(b));

cen_offset = 4; % Place fitting window this many bins after detected peak.
halfwin = 30;
%halfwin = ceil(halfwin);  % halfwin must be an integer.
%univar = (2*halfwin + 1)^2 / 12;  % Variance of a uniform distribution 
                                    % defined within fitting window.

thresholds_FA = b + 5*sqrt(b); % False alarm threshold
S = (pkvals > thresholds_FA) & (pkvals > threshold_absolute);

for ii = 1:size(dataCube, 1)
    cen = pkbins(ii)+cen_offset;
    if cen <= halfwin
        win_bins = (-cen+1):halfwin;
    elseif cen >= size(dataCube, 2)-halfwin
        win_bins = -halfwin:(size(dataCube, 2)-cen);
    else
        win_bins = -halfwin:halfwin;
    end
    signal = dataCube(ii, cen + win_bins);
    N = sum(signal);
    Delta = length(win_bins);
    bE = b(ii)*Delta;
    univar = Delta^2 / 12;
    E(ii) = max(N - bE, 0);
    a = max(E(ii)/N, 0);
    mu = sum(win_bins .* signal) ./ N;
    sigma2 = (1/(N-1))*sum(signal.*(win_bins - mu).^2);
    
    if a > 0
        p_mu(ii) = mu/a;
        p_w2(ii) = (1/a)*(sigma2 - (1-a)*univar) - (1-a)*p_mu(ii).^2;
        %p_w2(ii) = sigma2 - ((1-a)/a)*univar - (1-a)*p_mu(ii).^2;
        
        var_a = (b(ii)*Delta^2 / N^2) * (1/num_noise_bins + b(ii)/N);
        var_mu = 1/(12*N);
        
        var_p_mu(ii) = (1 / a^2) * (var_mu + (p_mu(ii)^2)*var_a);
        var_E(ii) = N + b(ii)*Delta^2/num_noise_bins;
        
        mupkcorr = interp1(mu_fit_LUT, pkbias_LUT, p_mu(ii));        
        p_pk(ii) = p_mu(ii) + mupkcorr + cen; 
        p_mu(ii) = p_mu(ii)+cen;
    end
end

