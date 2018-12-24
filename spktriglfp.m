function [avg_stlfp, sem_stlfp, accspk, stlfp_pow, stlfp_freq, band, ex] = spktriglfp(ex, varargin)
% spike triggered average lfp signal 
% [avg_stlfp, sem_stlfp, nspk] = spktriglfp( exSpk, exLFP)
%
% optional arguments are:
% 'time' - with following argument. restricting the lfp to the time before
%           and after a spike event. This also affects the number of spikes
%           because the surrounding LFP must be within the stimulus
%           presentation time. The default setting is +/-300ms.
% 'plot' - without following argument. Plots the spike triggered average
%           LFP
%
% @CL 09.04.2017


% warning(['I narrowed the analysis window to the window beginning 600ms after stimulus onset'...
%     'to 500ms before stimulus end in spktriglfp\getSpks. \n'...
%     'I do this to avoid dominant slow fluctuations in the beginning and the end for RC data. \n'])

% KK uses now data after 350ms...isn't it enough? (Chalk et al.,2010 used 256 - 512ms) 
% KK implemented firing rate normalization, as done by Chalk et al.,2010

p_flag = false;
wnd = 0.1; % was 0.3;  window before and after spike event to consider

%%% parse input
k = 1; 
while k<=length(varargin)
    switch varargin{k}
        case 'time'
            wnd = varargin{k+1};
        case 'plot'
            p_flag = true;  
    end
    k=k+1;
end

% preprocessing parameters
% nw = 2;
% nfft = 1024;
% Fs = 1000;
% [stlfp_pow, stlfp_freq] = pmtm(avg_stlfp, nw, nfft, Fs);
params = define_params;

%%% find spikes and estimate the lfp +/- time t around it
accstlfp = []; 
for t = 1:length(ex.Trials)
    [stlfp, nspk] = getstlfp4trial(ex.Trials(t), wnd);
    % trial-by-trial stLFP metrics
    ex.Trials(t).nspk = length(nspk);
    ex.Trials(t).mean_stLFP = nanmean(stlfp, 1);
    ex.Trials(t).sem_stLFP = nanstd(stlfp, [], 1)/sqrt(ex.Trials(t).nspk);
%     ex.Trials(t).mean_stLFP = nanmean(stlfp, 1)/ex.Trials(t).nspk;
%     ex.Trials(t).sem_stLFP = nanstd(stlfp/ex.Trials(t).nspk, [], 1)/sqrt(ex.Trials(t).nspk);
    [ex.Trials(t).stlfp_pow, ex.Trials(t).stlfp_freq] = mtspectrumc(ex.Trials(t).mean_stLFP, params);
    [ex.Trials(t).stlfp_delta, ex.Trials(t).stlfp_theta, ...
        ex.Trials(t).stlfp_alpha, ex.Trials(t).stlfp_beta, ex.Trials(t).stlfp_gamma]...
        = pow2band(ex.Trials(t).stlfp_freq, ex.Trials(t).stlfp_pow);
    
    % for trial-average
%     accstlfp = [accstlfp; ex.Trials(t).mean_stLFP];
        accstlfp = [accstlfp; stlfp];
end

%%% compute statistics
accspk = sum([ex.Trials.nspk]);
% avg_stlfp = nanmean(accstlfp, 1)/accspk; % average of the spike triggered lfp (stLFP)
avg_stlfp = nanmean(accstlfp, 1); % average of the spike triggered lfp (stLFP)
% avg_stlfp = avg_stlfp - mean(avg_stlfp);
sem_stlfp = nanstd(accstlfp, 0, 1)./sqrt(size(accstlfp,1)); % SEM of the stLFP

% % baseline correction
% avg_stlfp = avg_stlfp - mean(avg_stlfp);

% the PSD returned by pmtm is normalized per frequency unit
[stlfp_pow, stlfp_freq] = mtspectrumc(avg_stlfp, params);
[delta, theta, alpha, beta, gamma] = pow2band(stlfp_freq, stlfp_pow);
band = [delta, theta, alpha, beta, gamma];

%%% plot results
if p_flag
    
    [stimparam, vals] = getStimParam(ex);
    
    a1 = fill([-wnd:0.001:wnd, fliplr(-wnd:0.001:wnd)], ...
        [avg_stlfp - sem_stlfp , fliplr(avg_stlfp + sem_stlfp)], 'b');
    a1.FaceColor = [0.5 0.5 0.5]; a1.FaceAlpha = 0.4;
    a1.EdgeColor = 'w'; a1.EdgeAlpha = 0; hold on;
    
    plot(-wnd:0.001:wnd, avg_stlfp, ...
     'ButtonDownFcn', {@PlotAllSTA, -wnd:0.001:wnd, accstlfp}, ...
     'LineWidth', 2, ...
     'DisplayName', sprintf([stimparam '= %1.3f \n #spk: %1.0f'], vals, accspk),...
     'UserData', getFname(ex));
 
    xlabel('time rel:spike [s]');
    ylabel('avg LFP +/- SEM (\muV)');
    xlim([-wnd wnd]);
    legend('show');
    
end


%% Helper
function [stlfp, spikes] = getstlfp4trial(trial, wnd)
% lfp at spk +/- window wnd 

spikes = getSpks(trial, [0.35, 0]); % spikes within the stimulus presentation time
spikes = spikes{1};
stlfp = nan(length(spikes), length(-wnd:1/1000:wnd));

for i = 1:length(spikes)
    tspk = find(trial.LFP_prepro_time <= spikes(i), 1, 'last');
    tstrt = tspk-(wnd*1000);
    tend = tstrt+size(stlfp,2)-1;
    
    try
        % assign the spike centered LFP to the final matrix
        stlfp(i,:) = trial.LFP_prepro(tstrt:tend);
    catch
        disp('');        
    end
end

function [delta, theta, alpha, beta, gamma] = pow2band(freq, pow)
delta = nanmean(pow(freq > 0 & freq < 4));
theta = nanmean(pow(freq >= 4 & freq < 8));
alpha = nanmean(pow(freq >= 8 & freq <=13));
beta = nanmean(pow(freq > 13 & freq < 30));
gamma = nanmean(pow(freq >= 30 & freq <= 80));

%% Callback Function
function PlotAllSTA(source, ~,time, sta)

fprintf('plotting all spike triggered lfp \n')
nspikes = size(sta,1);

figure('Name', source.UserData);

for i = 1:nspikes
    plot(time, sta(i,:)); hold on;
end
set(findobj(gca, 'Type', 'Line'), 'Color', [0 0.2 0.2 0.2], 'LineWidth', 0.5);
plot(time, nanmean(sta), 'k'); hold on;
xlabel('time rel:spike [s]');
ylabel('\muV');

title(sprintf('single spike triggered LFP signals, #spikes %1d', nspikes));
crossl
