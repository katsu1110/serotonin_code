function [para, tuning] = LFPbyStm(ex)
% compute LFP parameters in the specified ex-file, split by types of
% stimuli
% - ex file must have LFP data, using 'loadCluster.m'
% - ex file must be filtered by 'filterLFP.m'

% stimulus type ==============
[stimparam, vals] = getStimParam(ex);
% blank_i =  vals > 1;
lenv = length(vals);

% initialization ==================
fs = 1000;
fsignal = 'MPsignal';

% structure
para.stm.param = stimparam;
para.stm.vals = vals;
para.ntr = length(ex.Trials);
para.bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
para.range = {[0.5,4], [4, 7], [8, 13], [14, 29], [30, 80]};
tuning.stmtype = stimparam;
tuning.res = [];
para.ts = ex.time;
para.f = ex.FREQ;
para.window = ex.period;
lenb = length(para.bands);
para.params = define_params;
para.wnd = 0.1; % was 0.3
% ampwindow = wndrange > -wnd &  wndrange  < wnd;
for i = 1:lenv    
    % trials with the unique stimulus
    trials = ex.Trials([ex.Trials.(stimparam)] == vals(i));
    lentr = length(trials);
    
    % analysis table
    % (stimulus, firing rate, LFP (uV))
    amat = nan(lentr, 3+lenb);
    amat(:,1) = vals(i);
    
    % firing rate per trial
    [~, spkc]  = getSpks(trials, [0 0]);
    amat(:,2) = spkc';
    
    % LFP trace averaged trials across the same stimulus
    lfpfull = vertcat(trials.(fsignal));
    para.lfp_stm.mean(i,:) = nanmean(lfpfull, 1);
    para.lfp_stm.sem(i,:) = nanstd(lfpfull, [], 1)/sqrt(lentr);
    amat(:,3) = nanmean(lfpfull(:, ex.time >= para.window{end}(1) ...
        & ex.time < para.window{end}(2)), 2);

    % spike-triggered average LFP & power frequency & Energy
    sta = [];
    para.spectrogram{i}.mean = zeros(length(para.f), length(para.ts));
    for n = 1:lentr
        % STA
        [spk, sc] = getSpks(trials(n), [para.window{end}(1) 0]);
        if sc > 0
            stlfp = getSTA(trials(n).(fsignal), trials(n).LFP_prepro_time, ...
                    spk{1}, para.wnd, fs);                
            sta = [sta; stlfp];
        end
        
        % power frequency
         for b = 1:lenb
             lfp = trials(n).(fsignal)(ex.time >= para.window{end}(1) ...
                & ex.time < para.window{end}(2));
             amat(n, 3+b) = bandpower(lfp, fs, para.range{b});
         end
        
         % energy
         para.spectrogram{i}.mean = para.spectrogram{i}.mean + trials(n).MPenergy;
    end
    if ~isempty(sta)
        para.stlfp.mean(i,:) = nanmean(sta, 1);
        para.stlfp.sd(i,:) = nanstd(sta, [], 1);
        para.stlfp.nspk(i) = size(sta, 1);
    else
        para.stlfp.mean(i,:) = nan(1, length(-para.wnd:1/fs:para.wnd));
        para.stlfp.sd(i,:) = para.stlfp.mean(i,:);
        para.stlfp.nspk(i) = 0;
    end
    %      % CL's correction
%        para.stlfp.avg_stlfp(i,:) = para.stlfp.avg_stlfp(i,:) ...
%            - mean(para.stlfp.avg_stlfp(i, wndrange < -0.06));

    tuning.res = [tuning.res; amat];
           
    % spectrogram
%     [para.period(u).spectrogram.S{i}, para.period(u).spectrogram.t{i}, para.period(u).spectrogram.f{i}] = ...
%         mtspecgramc(ph.period(u).lfpz', [0.1 0.01], para.params);
    para.spectrogram{i}.mean = para.spectrogram{i}.mean/lentr;

    % spike-LFP coherency
    spk = getSpks(trials, [para.window{end}(1) 0]);
    [para.coherence.C{i}, para.coherence.phi{i}, ...
        para.coherence.S12{i}, para.coherence.S1{i}, ...
        para.coherence.S2{i}, para.coherence.f{i}] = ...
            coherencycpt(lfpfull(:, ex.time >= para.window{end}(1) ...
                & ex.time < para.window{end}(2))', cell2struct(spk, 'spk', 1)', para.params);
end

%%
% tuning analysis
if length(para.stm.vals) > 2
    tuning.res(tuning.res(:,1) >= 1000, :) = [];
    for i = 2:size(tuning.res, 2)
        % information properties
        tuning.encoding{i-1} = encoding_tuning(tuning.res(:,1), tuning.res(:,i));
        % fit 
        tuning.fit{i-1} = fit_tuning(tuning.encoding{i-1}.mean, ...
            tuning.encoding{i-1}.std./sqrt(tuning.encoding{i-1}.ntr), ...
            tuning.encoding{i-1}.unistm, tuning.stmtype);
    end
else
    tuning.encoding = nan;
    tuning.fit = nan;
end