function para = stmLFP(ex, lfpfield, analysis)
%%
% LFP analysis by stimulus type
% INPUT: ex ... ex file after 'loadCluster.m'
%        lfpfield ... 'LFP_prepro', 'iLFP_prepro'
%        analysis ... 'all', 'sta'
%

% inputs =================
if nargin < 2; lfpfield = 'LFP_prepro'; end
if nargin < 3; analysis = {'all'}; end

% stimulus type ==============
[stimparam, vals] = getStimParam(ex);
% blank_i =  vals > 1;
lenv = length(vals);

% preset parameters ==================
fs = 2000;
para.stm.param = stimparam;
para.stm.vals = vals;
para.ts = ex.Trials(end).LFP_prepro_time;
para.wnd = 0.05; % was 0.3
if ex.exp.StimPerTrial == 4
    para.window = {[0.25 0.45]};
    stmdur = 0.45;
elseif ex.exp.StimPerTrial == 1
    para.window = {[0.5 2]};
    stmdur = 2;
end

% LFP as a function of stimulus
lfpfull = cell(1, lenv);
spk = cell(2, lenv);
spkc = cell(2, lenv);
lentr = zeros(1, lenv);
trials = cell(1, lenv);
para.spk_tu = cell(1, 2); 
for i = 1:lenv    
    % trials with the unique stimulus
    trials{i} = ex.Trials([ex.Trials.(stimparam)] == vals(i));
    lentr(i) = length(trials{i});
        
    % firing rate per trial in analysis window
    [spk{1, i}, spkc{1, i}]  = getSpks(trials{i}, [para.window{end}(1), ...
        stmdur - para.window{end}(2)]);    
    
    % firing rate per trial during stimulus presentation
    [spk{2, i}, spkc{2, i}]  = getSpks(trials{i}, [0 0]);
    
    % spikes elicited by the period (mean, SD, n)
    para.spk_tu{1}(i, :) = [mean(spkc{1, i}), std(spkc{1, i}), sum(spkc{1, i})];
    para.spk_tu{2}(i, :) = [mean(spkc{2, i}), std(spkc{2, i}), sum(spkc{2, i})];
    
    % removal of spike-related transient LFP
    for n = 1:lentr(i)
        lfpfull{i}(n, 1:length(trials{i}(n).(lfpfield))) = ...
            remove_spk(trials{i}(n).(lfpfield), ...
            trials{i}(n).LFP_prepro_time, spk{2, i}{n}, 0.003);
    end
    
%     % LFP trace averaged trials across the same stimulus
%     lfpfull{i} = vertcat(trials{i}.(lfpfield));
end
para.lfp = lfpfull;
para.ntr = lentr;

% Spike-triggered LFP =====================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'sta'))
%     ts = para.wnd*2;
%     movwin = [ts, ts/1000];
%     params = define_params(fs, ts);
    for i = 1:lenv
        sta = [];
        for n = 1:lentr(i)
            if spkc{1, i}(n) > 0                
                % compute a standard STA
                stlfp = getSTA(trials{i}(n).(lfpfield), trials{i}(n).LFP_prepro_time, ...
                    spk{1, i}{n}, para.wnd, fs);                
                sta = [sta; stlfp];
            end
        end
        if ~isempty(sta)
            % STA as a function of stimulus type
            para.sta.mean(i,:) = nanmean(sta, 1);
            para.sta.sd(i,:) = nanstd(sta, [], 1);
            para.sta.nspk(i) = size(sta, 1);
            
            % STA spectrogram
            [para.sta.s{i}, para.sta.f{i}, para.sta.t{i}, para.sta.p{i}] = ...
                spectrogram_frange(para.sta.mean(i,:), 95, fs, [0 100]);
%             [para.sta.S{i}, para.sta.t{i}, para.sta.f{i}] = ...
%                 mtspecgramc(sta', movwin, params);
        else
            para.sta.mean(i,:) = nan(1, length(-para.wnd:1/fs:para.wnd));
            para.sta.sd(i,:) = para.sta.mean(i,:);
            para.sta.nspk(i) = 0;
            para.sta.s{i} = nan;
            para.sta.f{i} = nan;
            para.sta.t{i} = nan;
            para.sta.p{i} = nan;
        end
    end   
end

% Spectrogram =============================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'spectrogram'))
    ts = para.window{end}(2) - para.window{end}(1);
    movwin = [ts, ts/100];
    params = define_params(fs, ts);
    for i = 1:lenv
%         [para.spectrogram.s{i}, para.spectrogram.f{i}, para.spectrogram.t{i}, para.spectrogram.p{i}]...
%                 = spectrogram_frange(mean(lfpfull{i}, 1), 90, fs, [0 100]);
        [para.spectrogram.S{i}, para.spectrogram.t{i}, para.spectrogram.f{i}] = ...
            mtspecgramc(lfpfull{i}', movwin, params);
    end   
end

% Spike-LFP coherency ========================
if ismember(1, contains(analysis, 'all')) || sum(contains(analysis, 'coherence'))
    wnd = [para.window{end}(1), para.window{end}(2), 0.003];
     params = define_params(fs, para.window{end}(2)-para.window{end}(1));
%     l = sum(ex.time >= para.window{end}(1) ...
%                 & ex.time <= para.window{end}(2));
    for i = 1:lenv
%         lfpmat = zeros(l, lentr(i));
%         spkmat = zeros(l, lentr(i));
%         for n = 1:lentr(i)
%             lfpmat(:, n) = lfpfull{i}(n, ex.time >= para.window{end}(1) ...
%                 & ex.time <= para.window{end}(2));
% %             spkmat(:, n) = b2c(spk{i}{n}, l, para.window{end});
%         end       
        [para.coherence.C{i}, para.coherence.phi{i}, ...
            para.coherence.S12{i}, para.coherence.S1{i}, ...
            para.coherence.S2{i}, para.coherence.f{i}] = ...
            spike_lfp_coherence(lfpfull{i}, trials{i}(end).LFP_prepro_time, spk{1, i}, wnd, params);
%         [para.coherence.cxy{i}, para.coherence.f{i}] = ...
%             mscohere(lfpmat, spkmat, hamming(round(l*0.2)), round(l*0.1), fs);
    end   
end

% % subfunction =========================================
function [C, phi, S12, S1, S2, f] = spike_lfp_coherence(lfpv, lfpt, spk, wnd, params)
% spike-LFP coherence analysis
ntr = length(spk);
lfp = lfpv(:, lfpt >= wnd(1) & lfpt <= wnd(2));
% lfpt = lfpt(lfpt >= wnd(1) & lfpt <= wnd(2)) - wnd(1);
s = struct('times', []);
for n = 1:ntr
    % adjust spike timings
    s(n).times = spk{n} - wnd(1);
    
%     % remove spike-component from LFP traces
%     lfp(n, :) = remove_spk(lfp(n, :), lfpt, spk{n}, wnd(3));
end
% get coherence by chronux function 
[C, phi, S12, S1, S2, f] = coherencycpt(lfp', s, params);    


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
nff = 256;

function [s, f, t, p] = spectrogram_frange(v, overlap, fs, frange)
[nsc, nov, nff] = stff_params(length(v), overlap);
[s, f, t, p] = spectrogram(v, nsc, nov, nff, fs);
outrange = f >= frange(1) & f <= frange(2);
s = s(outrange, :);
p = p(outrange, :);
f = f(outrange);
