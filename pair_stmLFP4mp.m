function para = pair_stmLFP4mp(ex0, ex2, exn0, exn2, lfpfield, thin, analysis)
%%
% LFP analysis by stimulus type with a pair of sessions
% INPUT: ex0, ex2 ... ex file after 'loadCluster.m'
%        exn0, exn2 ... exn file after 'MP_single.m'    
%        lfpfield ... 'LFP_prepro', 'iLFP_prepro'
%        thin ... 0, no thining, 1, thinning (firing rate correction)
%        analysis ... 'all', 'sta'
%

% inputs =================
if nargin < 5; lfpfield = 'LFP_prepro'; end
if nargin < 6; thin = 1; end
if nargin < 7; analysis = {'all'}; end

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
    para.window = {[0.2 0.45]};
    stmdur = 0.45;
    para.wnd = 0.05; % was 0.3
elseif ex0.exp.StimPerTrial == 1
    para.window = {[0.8 2]};
    stmdur = 2;
    para.wnd = 0.5; % was 0.3
end

% LFP & spikes as a function of stimulus
exs = {ex0, ex2};
exns = {exn0, exn2};
for d = 1:2
    % pupil vals
    for p = 1:length(exs{d}.Trials)
        % mean pupil size
        stpos = exs{d}.Trials(p).Eye_prepro.stpos;
        enpos = exs{d}.Trials(p).Eye_prepro.enpos;
        ps_temp = nanmean([exs{d}.Trials(p).Eye_prepro.psR(stpos:enpos); ...
            exs{d}.Trials(p).Eye_prepro.psL(stpos:enpos)], 1);

        % pupil size derivative
        dp_temp = nanmean([exs{d}.Trials(p).Eye_prepro.dpsR(stpos:enpos); ...
            exs{d}.Trials(p).Eye_prepro.dpsL(stpos:enpos)], 1);

        % assign
        if ex0.exp.StimPerTrial==1                
            exs{d}.Trials(p).psvals(1) = nanmean(ps_temp(end-round((enpos-stpos)/2)+1:end));               
            exs{d}.Trials(p).psvals(2) = max(dp_temp(end-round((enpos-stpos)/2)+1:end));
        else
            exs{d}.Trials(p).psvals(1) = nanmean(ps_temp);               
            exs{d}.Trials(p).psvals(2) = max(dp_temp);       
            if exs{d}.Trials(p).label_seq > 1
                for l = 1:exs{d}.Trials(p).label_seq-1
                    if exs{d}.Trials(p-l).label_seq < exs{d}.Trials(p-l+1).label_seq
                        exs{d}.Trials(p-l).psvals = exs{d}.Trials(p).psvals;
                    end
                end
            end
        end
    end
    
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
        mptrs = exns{d}.Trials([exs{d}.Trials.(stimparam)] == vals(i));
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
        
        % swap preprocessed LFP
        for n = 1:para.cond(d).ntr(i)
            % replace with MPs
            para.cond(d).trials{i}(n).(lfpfield) = mptrs(n).signal;
            para.cond(d).trials{i}(n).energy = mptrs(n).energy;
            para.cond(d).lfpfull{i}(n, 1:length(para.cond(d).trials{i}(n).(lfpfield))) = ...
                para.cond(d).trials{i}(n).(lfpfield);
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

% Synchronized index ===============================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'si'))
    for d = 1:2
        for i = 1:lenv
            for n = 1:para.cond(d).ntr(i)
                para.cond(d).si{i}(n) = spktrain2psi(para.cond(d).spk{j, i}{n}, para.window{end}, 50, fs);
            end
        end
    end       
end

% Spike-triggered LFP =====================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'sta'))
    for d = 1:2
        for i = 1:lenv
            sta = [];
            for n = 1:para.cond(d).ntr(i)
                try      
                    % compute a standard STA
                    stlfp = getSTA(para.cond(d).trials{i}(n).(lfpfield), ...
                        para.cond(d).trials{i}(n).LFP_prepro_time, ...
                        para.cond(d).spk{j, i}{n}, para.wnd, fs);                
                    sta = [sta; stlfp];
                catch
                    continue
                end
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

% Tuning ================================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'tuning'))
    if ~strcmp(para.stm.param, 'rc')
        % unique params
        bandnames = {'spikes', 'theta', 'alpha', 'beta', 'slow_gamma', 'fast_gamma'};
        bandrange = {[0, 7], [8, 13], [14, 24], [25 40], [45, 70]};
        lenb = length(bandrange);
                
        % data extraction
        for d = 1:2
            tumat = nan(sum(para.cond(d).ntr), 2+lenb);
            c = 1;
            for i = 1:lenv
                if para.stm.vals(i) < 1000
                    % stimulus
                    tumat(c:c+para.cond(d).ntr(i)-1, 1) = para.stm.vals(i);

                    % spikes 
                    tumat(c:c+para.cond(d).ntr(i)-1, 2) = para.cond(d).spkc{1, i};

                    % spectra
                    f_res = exns{1}.freq;
                    t = exns{1}.MPtime;
                    for m = 1:para.cond(d).ntr(i)
                        for b = 1:lenb
                            S_res = para.cond(d).trials{i}(m).energy(...
                                f_res >= bandrange{b}(1) & f_res <= bandrange{b}(2), ...
                                para.window{end}(1)<=t & t<=para.window{end}(2));
                             tumat(c:c+m-1, 2+b) = 10*log10(nanmean(S_res(:)));
                        end
                    end
                end
                
                c = c + para.cond(d).ntr(i);
            end   

            % encoding analysis
            tumat(isnan(tumat(:, 1)), :) = [];
            for b = 1:lenb+1
                try
                    para.cond(d).tuning.(bandnames{b}) = encoding_tuning(...
                        tumat(:, 1), tumat(:, 1+b), para.stm.param);  
                catch
                    disp(['Tuning analysis for ' bandnames{b} ' skipped.'])
                    para.cond(d).tuning.(bandnames{b}) = nan;
                end
            end
        end
    end
end


% Spectrogram =============================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'spectrogram'))
    for d = 1:2
        for i = 1:lenv
            para.cond(d).spectrogram.f{i} = exns{1}.freq;
            para.cond(d).spectrogram.t{i} = exns{1}.MPtime;
            for m = 1:para.cond(d).ntr(i)
                if m==1
                    para.cond(d).spectrogram.S{i} = para.cond(d).trials{i}(m).energy;
                else
                    if sum(isnan(para.cond(d).trials{i}(m).energy(:)))==0
                        para.cond(d).spectrogram.S{i} = para.cond(d).spectrogram.S{i} + para.cond(d).trials{i}(m).energy;
                    end
                end
            end
            para.cond(d).spectrogram.S{i} = para.cond(d).spectrogram.S{i}/para.cond(d).ntr(i);
        end   
    end
end

% Spike-LFP coherency ========================
if ismember(1, contains(analysis, 'all')) || sum(contains(analysis, 'coherence'))
    ts = 0.2;
    params = define_params(fs, ts, 1);
    for d = 1:2
        for i = 1:lenv
            % coherence analysis by chronux toolbox
            [para.cond(d).coherence.C{i}, para.cond(d).coherence.phi{i}, ...
                para.cond(d).coherence.S12{i}, para.cond(d).coherence.S1{i}, ...
                para.cond(d).coherence.S2{i}, para.cond(d).coherence.f{i}] = ...
                spike_lfp_coherence(para.cond(d).lfpfull{i}, para.cond(d).trials{i}(end).LFP_prepro_time, ...
                para.cond(d).spk{j, i}, para.window{end}, params);
            
%             % trial average, except for phase
%             para.cond(d).coherence.C{i} = mean(para.cond(d).coherence.C{i}, 2);
%             para.cond(d).coherence.S12{i} = mean(para.cond(d).coherence.S12{i}, 2);
%             para.cond(d).coherence.S1{i} = mean(para.cond(d).coherence.S1{i}, 2);
%             para.cond(d).coherence.S2{i} = mean(para.cond(d).coherence.S2{i}, 2);
%             para.cond(d).coherence.f{i} = mean(para.cond(d).coherence.f{i}, 1);
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
