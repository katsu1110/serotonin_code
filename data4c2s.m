function [para, stlfp0, stlfp1, data] = data4c2s(ex0, ex2, lfpfield, thin, z)
%%
% format LFP data from the RC experiments such that they can be used in
% 'c2s' analysis
%
% source code for c2s : 'https://github.com/lucastheis/c2s'
%
% INPUT: ex0, ex2 ... ex file after 'loadCluster.m'
%        lfpfield ... 'LFP_prepro', 'iLFP_prepro'
%        thin ... 0, no thining, 1, thinning (firing rate correction)
%        analysis ... 'all', 'sta'
%        z ... z-score (0 or 1)
%
% OUTPUT: para ... struct; session information
%         stlfp1 ... spike-triggered LFPs
%         stlfp0 ... LFP segments without spikes in the center
%         data ... cell array; something as follows to match the convention of 'c2s'
% data{1} =
%      calcium: [1x71985 double] .... in our case, preprocessed LFP
%       spikes: [1x71985 uint16] .... spikes
%          fps: 99.9998          .... sampling rate 
%     cell_num: 1                .... 1, baseline; 2, drug
%

% inputs =================
if nargin < 3; lfpfield = 'LFP_prepro'; end
if nargin < 4; thin = 1; end
if nargin < 5; z = 1; end

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
    para.wnd = 0.07; % was 0.3
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
if z==1
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
    me = zeros(1, 1); sd = ones(1, 1);
end
for d = 1:2
    lenlfp = length(para.cond(d).trials{end}(end).(lfpfield));
    for i = 1:lenv    
        for n = 1:para.cond(d).ntr(i)
            para.cond(d).trials{i}(n).(lfpfield) = (...
                para.cond(d).trials{i}(n).(lfpfield) - me(1))/sd(1);
            para.cond(d).lfpfull{i}(n, 1:lenlfp) = para.cond(d).trials{i}(n).(lfpfield);
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

% generate 'data' for c2s analysis
data = cell(1, 2);
for d = 1:2
    % sampling rate
    data{d}.fps = fs;
    
    % condition
    data{d}.cell_num = d;
    
    % initialize
    data{d}.spike_times = [];
    data{d}.calcium = [];
    trt = 0;
    for i = 1:lenv 
        for n = 1:para.cond(d).ntr(i)
            % spike times (ms)
            spkt = 1000*(para.cond(d).spk{3, i}{n} - para.window{1}(1)) + trt;
            data{d}.spike_times = [data{d}.spike_times, spkt'];
            trt = trt + 1000*(para.window{1}(2) - para.window{1}(1));

            % LFP
            lfp = para.cond(d).lfpfull{i}(n, para.ts >= para.window{1}(1)...
                & para.ts <= para.window{1}(2));
            data{d}.calcium = [data{d}.calcium, lfp];
        end
    end
end

% stLFPs
stlfp0 = cell(1,2);
stlfp1 = cell(1,2);
ts = para.window{1}(1):0.0001:para.window{1}(2);
for d = 1:2
    for i = 1:lenv
        for n = 1:para.cond(d).ntr(i)
            % compute a standard STA
            stlfp = getSTA(para.cond(d).trials{i}(n).(lfpfield), ...
                para.cond(d).trials{i}(n).LFP_prepro_time, ...
                para.cond(d).spk{j, i}{n}, para.wnd, fs);
            
            % LFP without spikes
            t_temp = ts';
            tout = [];
            for k = 1:length(para.cond(d).spk{j, i}{n})
                tout = [tout, para.cond(d).spk{j, i}{n}(k)-...
                    para.wnd:0.0001:para.cond(d).spk{j, i}{n}(k)+para.wnd];
            end
            t_temp(ismember(t_temp, tout)) = [];
            stlfp_null = getSTA(para.cond(d).trials{i}(n).(lfpfield), ...
                para.cond(d).trials{i}(n).LFP_prepro_time, ...
                t_temp, para.wnd, fs);
            
            if ~isnan(stlfp)
                % stack stLFP
                stlfp1{d} = [stlfp1{d}; stlfp]; 
                stlfp0{d} = [stlfp0{d}; stlfp_null];
            end
        end            
    end  
end     

% remove irrelevant fields
rmf = {'ts', 'window', 'cond'};
for r = 1:length(rmf)
    para = rmfield(para, rmf{r});
end

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