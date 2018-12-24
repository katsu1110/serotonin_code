function [prshat, ff] = fit_flexiblemodulatedpoiss(spikec, mdltyp)
%%
% fit the flexible modulated poisson model to spike count data
%
% INPUT: spikec ... spike counts: trials x discrete stimuli
%              mdltyp ... model type: 'exp', 'softrect', or 'softrectp'
%
% OUTPUT: prshat ... (1) stimulus dependent term, (2) noise variance, (3) power exponent (only for 'softrectpow')
%                 ff ... fano factor
%
% NOTE that 'https://github.com/adamshch/flexibleModulatedPoisson' must be
% cloned beforehand on your path
%

if nargin < 2; mdltyp = 'exp'; end
n_ori = size(spikec, 2); % the number of discrete stimuli
opts   = optimset('largescale','off','maxfunevals',1e5);  % Set optimization paraemters for fminunc

switch mdltyp
    case 'softrect'
        prs0 = zeros(n_ori+1, 1);
        % Define negative log-likelihood (for a soft-rectification nonlinearity rased to a power)
        negLsrp = @(x) negLfun_latentPoiss(x,@logli_latentPoiss_softrect, spikec);
    case 'exp'
        prs0 = zeros(n_ori+1, 1);
        % Define negative log-likelihood (for an exponential nonlinearity rased to a power)
        negLsrp = @(x) negLfun_latentPoiss(x,@logli_latentPoiss_exp, spikec);
    case 'softrectp'
        prs0 = zeros(n_ori+2, 1);
        % Define negative log-likelihood (for a soft-rectification power nonlinearity rased to a power)
        negLsrp = @(x) negLfun_latentPoiss(x,@logli_latentPoiss_softrectpow, spikec);
end
% Run fminunc to find the parameters (for a soft-rectification nonlinearity rased to a power)
prshat = fminunc(negLsrp,prs0,opts);  

% fano factor
ff = nan(1, n_ori);
for n = 1:n_ori 
    switch mdltyp
        case 'softrectp'
            % flexible power exponent. Delta method.
%             error('not sure how to compute mean and variance here.')
            me = exp(prshat(n) + 0.5*prshat(end));
            va = me + exp(prshat(end) - 1)*me^2;
        otherwise
            % exponential nonlinearity
            me = exp(prshat(n) + 0.5*prshat(end));
            va = me + (exp(exp(prshat(end))) - 1)*(me^2);
%             va = me + prshat(end)*(exp(prshat(n)))^2;
    end
    ff(n) = va/me;
end