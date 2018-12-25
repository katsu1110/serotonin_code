function para = pair_stmLFP(ex0, ex2, lfpfield, thin, analysis, z)
%%
% LFP analysis by stimulus type with a pair of sessions
% INPUT: ex0, ex2 ... ex file after 'loadCluster.m'
%        lfpfield ... 'LFP_prepro', 'iLFP_prepro'
%        thin ... 0, no thining, 1, thinning (firing rate correction)
%        analysis ... 'all', 'sta'
%        z ... z-score (0 or 1)
%

% inputs =================
if nargin < 3; lfpfield = 'LFP_prepro'; end
if nargin < 4; thin = 1; end
if nargin < 5; analysis = {'all'}; end
if nargin < 6; z = 1; end

% stimulus type ==============
[stimparam, vals] = getStimParam(ex0);
lenv = length(vals);
if lenv==1 && strcmp(stimparam, 'or')
    para.stm.param = 'rc';
else
    para.stm.param = stimparam;
end
para.stm.vals = vals;

% drug & animal ==================
try
    para.ismango = contains(ex2.Header.fileName, 'Mango');
    para.is5ht = contains(ex2.Header.fileName, '5HT');
catch
    para.ismango = nan;
    para.is5ht = nan;
end
    
% preset parameters ==================
fs = 1000;
para.ts = ex0.Trials(end).LFP_prepro_time;
if ex0.exp.StimPerTrial == 4
    para.window = {[0.25 0.45]};
    stmdur = 0.45;
    para.wnd = 0.05; % was 0.3
elseif ex0.exp.StimPerTrial == 1
    para.window = {[0.8 2]};
    stmdur = 2;
    para.wnd = 0.5; % was 0.3
end

% LFP & spikes as a function of stimulus
exs = {ex0, ex2};
for d = 1:2
    % initialization
    para.cond(d).lfpfull = cell(1, lenv);
    para.cond(d).spk = cell(3, lenv);
    para.cond(d).spkc = cell(3, lenv);
    para.cond(d).ntr = zeros(1, lenv);
    para.cond(d).trials = cell(1, lenv);
    para.cond(d).spk_tu = cell(1, 3); 
    
    % stimulus types
    for i = 1:lenv    
        % trials with the unique stimulus
        para.cond(d).trials{i} = exs{d}.Trials([exs{d}.Trials.(stimparam)] == vals(i));
        para.cond(d).ntr(i) = length(para.cond(d).trials{i});        

        % firing rate per trial in analysis window
        [para.cond(d).spk{1, i}, para.cond(d).spkc{1, i}]  = ...
            getSpks(para.cond(d).trials{i}, [para.window{end}(1), ...
            stmdur - para.window{end}(2)]);    

        % firing rate per trial during stimulus presentation
        [para.cond(d).spk{2, i}, para.cond(d).spkc{2, i}]  = getSpks(para.cond(d).trials{i}, [0 0]);

        % elicited spikes (mean, SD, n)
        dur = para.window{end}(2) - para.window{end}(1);
        para.cond(d).spk_tu{1}(i, :) = [mean(para.cond(d).spkc{1, i})/dur, ...
            std(para.cond(d).spkc{1, i})/dur, sum(para.cond(d).spkc{1, i})];        
        para.cond(d).spk_tu{2}(i, :) = [mean(para.cond(d).spkc{2, i})/stmdur, ...
            std(para.cond(d).spkc{2, i})/stmdur, sum(para.cond(d).spkc{2, i})];    
        
        % preprocessed LFP
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).lfpfull{i}(n, 1:length(para.cond(d).trials{i}(n).(lfpfield))) = ...
                para.cond(d).trials{i}(n).(lfpfield);
        end
    end
end

% zscoring
if z
    lfpall = [];
    for d = 1:2
        for i = 1:lenv    
            for n = 1:para.cond(d).ntr(i)
                lfpall = [lfpall, para.cond(d).trials{i}(n).(lfpfield)];
            end
        end
    end
    me = nanmean(lfpall);
    sd = nanstd(lfpall);
else
    me = 0; sd = 1;
end
for d = 1:2
    lenlfp = length(para.cond(d).trials{end}(end).(lfpfield));
    for i = 1:lenv    
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).lfpfull{i}(n, 1:lenlfp) = (...
                para.cond(d).trials{i}(n).(lfpfield) - me)/sd;
        end
    end
end
    
% 'thinning' to control the difference in firing rate 
j = 1;
if thin==1
    j = 3;
    a = [1, 2];
    for i = 1:lenv
        % rate ratio
        alpha = para.cond(2).spk_tu{1}(i, 1)/para.cond(1).spk_tu{1}(i, 1);
        if alpha == 1
            para.cond(1).spk_tu{3}(i, :) = para.cond(1).spk_tu{1}(i, :);
            para.cond(2).spk_tu{3}(i, :) = para.cond(2).spk_tu{1}(i, :);
            para.cond(1).spk{3, i} = para.cond(1).spk{1, i};
            para.cond(1).spkc{3, i} = para.cond(1).spkc{1, i};
            para.cond(2).spk{3, i} = para.cond(2).spk{1, i};
            para.cond(2).spkc{3, i} = para.cond(2).spkc{1, i};
            continue
        elseif alpha < 1
            d = 1;
        elseif alpha > 1
            d = 2;
            alpha = 1/alpha;
        end
        para.cond(d).spk{3, i} = para.cond(d).spk{1, i};
        para.cond(d).spkc{3, i} = para.cond(d).spkc{1, i};
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).spk{3, i}{n} = thinning(para.cond(d).spk{1, i}{n}, alpha);
            para.cond(d).spkc{3, i}(n) = length(para.cond(d).spk{3, i}{n});
        end
        b = a(~ismember(a, d));
        para.cond(b).spk{3, i} = para.cond(b).spk{1, i};
        para.cond(b).spkc{3, i} = para.cond(b).spkc{1, i};

        % elicited spikes (mean, SD, n) after thinning
        para.cond(d).spk_tu{3}(i, :) = [mean(para.cond(d).spkc{3, i})/dur, ...
                std(para.cond(d).spkc{3, i})/dur, sum(para.cond(d).spkc{3, i})];     
        para.cond(b).spk_tu{3}(i, :) = para.cond(b).spk_tu{1}(i, :);     
    end
end

% trial-by-trial analysis ====================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'mat'))
    collab = {'spike count (su)', 'spike count (mu)', 'mean pupil size', ...
                'derivative pupil size', 'lfp response','low-freq (<10Hz)', 'gamma-freq (>40Hz)'};
    lenl = length(collab);
    ts = 0.3; % 3 tapers
    params = define_params(fs, ts, 0);
    for d = 1:2
        para.cond(d).collab = collab; 
        for i = 1:lenv            
            para.cond(d).mat{i} = nan(para.cond(d).ntr(i), lenl);
            
            % MU spike counts
            [~, spc0]  = getSpks(para.cond(d).trials{i}, [0 0], 'oSpikes');
            
            % LFP power
            [S, f] = mtspectrumc(para.cond(d).lfpfull{i}(:, ...
                para.cond(d).trials{i}(end).LFP_prepro_time > para.window{end}(1) ...
                & para.cond(d).trials{i}(end).LFP_prepro_time <= para.window{end}(2)), ...
                params);
            for n = 1:para.cond(d).ntr(i)
                % spike count (su)
                para.cond(d).mat{i}(n, 1) = para.cond(d).spkc{2, i}(n);
                
                % spike count (mu)
                para.cond(d).mat{i}(n, 2) = spc0(n);
                
                % mean pupil size
                stpos = para.cond(d).trials{i}(n).Eye_prepro.stpos;
                enpos = para.cond(d).trials{i}(n).Eye_prepro.enpos;
                ps_temp = nanmean([para.cond(d).trials{i}(n).Eye_prepro.psR; ...
                    para.cond(d).trials{i}(n).Eye_prepro.psL], 1);
                para.cond(d).mat{i}(n, 3) = nanmean(ps_temp(end-round((enpos-stpos)/4)+1:end));
                
                % pupil size derivative
                ps_temp = nanmean([para.cond(d).trials{i}(n).Eye_prepro.dpsR; ...
                    para.cond(d).trials{i}(n).Eye_prepro.dpsL], 1);
                para.cond(d).mat{i}(n, 4) = nanmean(ps_temp(end-round((enpos-stpos)/4)+1:end));
                
                % lfp response
                lfpseg = para.cond(d).trials{i}(n).(lfpfield)(para.cond(d).trials{i}(n).LFP_prepro_time > 0 ...
                    & para.cond(d).trials{i}(n).LFP_prepro_time <= 0.2);
                para.cond(d).mat{i}(n, 5) = max(lfpseg) - min(lfpseg);
                
                % low-freq power
                para.cond(d).mat{i}(n, 6) = nanmean(10*log10(S(f <= 10, n)), 1);
                
                % gamma power
                para.cond(d).mat{i}(n, 7) = nanmean(10*log10(S(f >= 40, n)), 1);
            end
        end  
    end     
end

% Spike-triggered LFP ====================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'sta'))
    for d = 1:2
        for i = 1:lenv
            sta = [];
            for n = 1:para.cond(d).ntr(i)
                % compute a standard STA
                stlfp = getSTA(para.cond(d).trials{i}(n).(lfpfield), ...
                    para.cond(d).trials{i}(n).LFP_prepro_time, ...
                    para.cond(d).spk{j, i}{n}, para.wnd, fs);                
                sta = [sta; stlfp];
            end
            if ~isempty(sta)
                % STA as a function of stimulus type
                para.cond(d).sta.mean(i,:) = nanmean(sta, 1);
                para.cond(d).sta.sd(i,:) = nanstd(sta, [], 1);
                para.cond(d).sta.nspk(i) = size(sta, 1);

                % STA spectrogram
                [para.cond(d).sta.s{i}, para.cond(d).sta.f{i}, para.cond(d).sta.t{i}, para.cond(d).sta.p{i}] = ...
                    spectrogram_frange(para.cond(d).sta.mean(i,:), 95, fs, [0 100]);
            else
                para.cond(d).sta.mean(i,:) = nan(1, length(-para.wnd:1/fs:para.wnd));
                para.cond(d).sta.sd(i,:) = para.cond(d).sta.mean(i,:);
                para.cond(d).sta.nspk(i) = 0;
                para.cond(d).sta.s{i} = nan;
                para.cond(d).sta.f{i} = nan;
                para.cond(d).sta.t{i} = nan;
                para.cond(d).sta.p{i} = nan;
            end
        end  
    end     
end


% % Tuning ================================
% if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'tuning'))
%     if ~strcmp(para.stm.param, 'rc')
%         % unique params
%         bandnames = {'spikes', 'theta', 'alpha', 'beta', 'slow_gamma', 'fast_gamma'};
%         bandrange = {[0, 7], [8, 13], [14, 24], [25 40], [45, 70]};
%         lenb = length(bandrange);
%                 
%         % data extraction
%         for d = 1:2
%             tumat = nan(sum(para.cond(d).ntr), 2+lenb);
%             c = 1;
%             for i = 1:lenv
%                 if para.stm.vals(i) < 1000
%                     % stimulus
%                     tumat(c:c+para.cond(d).ntr(i)-1, 1) = para.stm.vals(i);
% 
%                     % spikes 
%                     tumat(c:c+para.cond(d).ntr(i)-1, 2) = para.cond(d).spkc{1, i};
% 
%                     % spectra
%                     f_res = exns{1}.freq;
%                     t = exns{1}.MPtime;
%                     for m = 1:para.cond(d).ntr(i)
%                         for b = 1:lenb
%                             S_res = para.cond(d).trials{i}(m).energy(...
%                                 f_res >= bandrange{b}(1) & f_res <= bandrange{b}(2), ...
%                                 para.window{end}(1)<=t & t<=para.window{end}(2));
%                              tumat(c:c+m-1, 2+b) = 10*log10(nanmean(S_res(:)));
%                         end
%                     end
%                 end
%                 
%                 c = c + para.cond(d).ntr(i);
%             end   
% 
%             % encoding analysis
%             tumat(isnan(tumat(:, 1)), :) = [];
%             for b = 1:lenb+1
%                 try
%                     para.cond(d).tuning.(bandnames{b}) = encoding_tuning(...
%                         tumat(:, 1), tumat(:, 1+b), para.stm.param);  
%                 catch
%                     disp(['Tuning analysis for ' bandnames{b} ' skipped.'])
%                     para.cond(d).tuning.(bandnames{b}) = nan;
%                 end
%             end
%         end
%     end
% end


% Spectrogram =============================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'spectrogram'))
    ts = 0.3; % 3 tapers
    params = define_params(fs, ts, 1);
    movwin = [ts, 0.01];
    for d = 1:2
        for i = 1:lenv
            % spectrogram by chronux toolbox
            [para.cond(d).spectrogram.S{i}, para.cond(d).spectrogram.t{i}, ...
                para.cond(d).spectrogram.f{i}] = mtspecgramc(para.cond(d).lfpfull{i}', movwin, params);
        end   
    end
end

% Spike-LFP coherency ========================
if ismember(1, contains(analysis, 'all')) || sum(contains(analysis, 'coherence'))
    ts = 0.3; % 3 tapers
    params = define_params(fs, ts, 0);
    for d = 1:2
        for i = 1:lenv
            % coherence analysis by chronux toolbox
            [~, ~, para.cond(d).coherence.S12{i}, para.cond(d).coherence.S1{i}, ...
                para.cond(d).coherence.S2{i}, para.cond(d).coherence.f{i}] = ...
                spike_lfp_coherence(para.cond(d).lfpfull{i}, para.cond(d).trials{i}(end).LFP_prepro_time, ...
                para.cond(d).spk{j, i}, para.window{end}, params);
            
            % trial average to deal with nans
            para.cond(d).coherence.S12{i} = squeeze(nanmean(para.cond(d).coherence.S12{i}, 2));
            para.cond(d).coherence.S1{i} = squeeze(nanmean(para.cond(d).coherence.S1{i}, 2));
            para.cond(d).coherence.S2{i} = squeeze(nanmean(para.cond(d).coherence.S2{i}, 2));
            para.cond(d).coherence.f{i} = nanmean(para.cond(d).coherence.f{i}, 1)';
            
            % trial average coherence and phase
            C12 = para.cond(d).coherence.S12{i}./sqrt(para.cond(d).coherence.S1{i}.*para.cond(d).coherence.S2{i});
            para.cond(d).coherence.C{i} = abs(C12);
            para.cond(d).coherence.phi{i} = angle(C12);
        end   
    end
end

% remove fields for memory =====================
rms = {'spk', 'spkc', 'trials'};
for i = 1:length(rms)
    para.cond = rmfield(para.cond, rms{i});
end

% % subfunction =========================================
function [C, phi, S12, S1, S2, f] = spike_lfp_coherence(lfpv, lfpt, spk, wnd, params)
% spike-LFP coherence analysis
ntr = length(spk);
% % binary vector for spikes
% lent = length(lfpt);
% s = struct('times', []);
% s = zeros(lent, ntr);
% for n = 1:ntr
%     s(n).times = zeros(lent, 1);
%     nspk = length(spk{n});
%     for k = 1:nspk
%         [~, idx] = min(abs(lfpt - spk{n}(k)));
% %         s(n).times(idx) = 1;
%         s(idx, n) = 1;
%     end
%     s(n).times = s(n).times(lfpt >= wnd(1) & lfpt <= wnd(2));
% end
% s = s(lfpt >= wnd(1) & lfpt <= wnd(2), :);
lfp = lfpv(:, lfpt >= wnd(1) & lfpt <= wnd(2));
s = struct('times', []);
for n = 1:ntr
    % adjust spike timings
    s(n).times = spk{n} - wnd(1);
end
% get coherence by chronux function 
% [C, phi, S12, S1, S2, f] = coherencyc(lfp', s, params);    
[C, phi, S12, S1, S2, f] = coherencycpt(lfp', s, params, 1); 

function spk = thinning(spk, fr_ratio)
%%
% perform thinning of spikes
% INPUT: spk ... spike times (1D)
%              fr_ratio ... ratio of mean firing rate across two conditions
%                              (u_small/u_large)
%              wnd ... window to apply thinning
%

if isempty(spk)
    return
end

r = rand(1, length(spk));
spk(r <= 1 - fr_ratio) = []; 

% function sc = b2c(spk, l, wnd)
% sc = zeros(l, 1);
% t = linspace(wnd(1), wnd(2), l);
% for i = 1:length(spk)
%     [~, idx] = min(abs(t - spk(i)));
%     sc(idx) = 1;
% end

function [nsc, nov, nff] = stff_params(L, overlap)
nsc = floor(L/18);
nov = floor(nsc/(100/overlap));
% nff = max(256, 2^nextpow2(nsc));
nff = 1024;

function [s, f, t, p] = spectrogram_frange(v, overlap, fs, frange)
[nsc, nov, nff] = stff_params(length(v), overlap);
[s, f, t, p] = spectrogram(v, nsc, nov, nff, fs);
outrange = f >= frange(1) & f <= frange(2);
s = s(outrange, :);
p = p(outrange, :);
f = f(outrange);
