function k = coherence_rate_adjusted(fr1, fr2, sps, hpn, dt)
%%
% yields the rate-adjustment factor (Aoi et al., 2015) for Spike-LFP
% coherence
% INPUT: fr1, fr2 ... mean firing rate in the two conditions (fr1 > fr2)
%        sps ... spike power spectrum in the condition 1
%        hpn ... homogeneous poisson noise
%               (intercept with the affine transformation)
%        dt ... size of time step
%

if fr2 > fr1
    error('fr1 must be > fr2.')
end
if nargin < 3
    error('At least mean firing rates in two conditions and spike power spectrum must be used as input arguments.'); 
end
if nargin < 4; hpn = 0; end
if nargin < 5; dt = 1; end

alpha = fr2/fr1;
adjusted = ((dt)^2)*((1/alpha - 1)*fr1 + hpn/alpha^2);
k = 1/sqrt(1 + adjusted/sps);
