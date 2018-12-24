function prshat_srp = fit_soft_rect_p(ex, spkfield, wnd, detrendbin)
%% 
% fit 'soft_rect_p' by Charles et al., 2018; 'Dethroning the Fano Factor: A Flexible, Model-Based
% Approach to Partitioning Neural Variability'
%
% INPUT: ex ... ex-file
%              spkfield ... 'Spikes' or 'oSpikes'                   
%              wnd ... analysis window ([0.8 0] in default)
%              detremdbin ... perform local detrending across trials with
%              specfied bin (trials) size
%
% OUTPUT: prshat_srp .... fitting results
% prshat_srp(1:end-2);    % Extract optimized mean parameter (power soft-rectification)
% exp(prshat_srp(end-1)); % Extract optimized variance parameter (power soft-rectification)
% exp(prshat_srp(end));   % Extract optimized variance parameter (power soft-rectification)
%
%

%%
% trial-by-trial firing rate
ex.Trials = ex.Trials(abs([ex.Trials.Reward]) > 0);
[~, spkc] = getSpks(ex.Trials, wnd, spkfield);

%%
% detrending
if detrendbin > 0
    % TODO
    
end

%%
% fit 'soft_rect_p'
% size(spkc)
negLsrp = @(x) negLfun_latentPoiss(x,@logli_latentPoiss_softrectpow, spkc');% Define negative log-likelihood (for a soft-rectification nonlinearity rased to a power)
n_ori = 1; % assume only one orientation for the RC experiment
opts   = optimset('display','iter','largescale','off','maxfunevals',1e5);  % Set optimization paraemters for fminunc
prs0   = zeros(n_ori+2,1);  % Initialize log-likelihood (for an exponential/softrect nonlinearity this is the number of orientations plus 1 variance parameter and 1 power parameter)
prshat_srp = fminunc(negLsrp,prs0,opts);  % Run fminunc to find the parameters (for a soft-rectification nonlinearity rased to a power)
