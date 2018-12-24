function plot_LFPs(serdata, stmtype)
%% plot LFP data across sessions
% INPUT: serdata
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_drug.mat')
%
% ++++++++++++++++++++++++++++++++

close all;
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc'))

%% 
% pair type
pairtype = serdata.session(459).results.pairtype;
switch pairtype
    case 'drug'
        names = {'baseline', 'drug'};
    case 'sc'
        names = {'low sc', 'high sc'};
    case 'ps'
        names = {'small ps', 'large ps'};
    case 'ps_drug'
        names = {'small ps drug', 'large ps drug'};
end

% stimulus type
row = [];
is5ht = [];
for i = 1:length(serdata.session)    
    if serdata.session(i).exist==1 && serdata.session(i).goodunit==1 ...
        && isstruct(serdata.session(i).results) && strcmp(stmtype, serdata.session(i).results.stimulus)  
        row = [row, i];        
        if strcmp(serdata.session(i).results.drugname, '5HT')
            is5ht = [is5ht, 1];
        else
            is5ht = [is5ht, 0];
        end
    end
end

lenr = length(row);
disp(['Analized pairs of sessions: ' num2str(lenr)])
disp(['5HT sessions: ' num2str(sum(is5ht==1))])
disp(['NaCl sessions: ' num2str(sum(is5ht==0))])

serdata.session = serdata.session(row);

%%
% data extraction
para.is5ht = is5ht;
para.pairtype = pairtype;
para.stmtype = stmtype;
para.names = names;
if strcmp(stmtype, 'rc')
    len_trace_lfp = 2000;
else
    len_trace_lfp = 450;
end
%     len_trace_lfp = size(serdata.session(row(1)).results.(fieldname).cond(1).lfpstm.lfp_stm.mean, 2);
%     len_trace_lfp = 2000;
xfreq = serdata.session(1).results.cond(1).lfpstm.f;
len_pow_lfp = length(xfreq);
for k = 1:2
    para.cond(k).lfpstm.lfp.trace = nan(lenr, len_trace_lfp);
    para.cond(k).lfpstm.lfp.power = nan(lenr, len_pow_lfp);
end
for i = 1:lenr
    for k = 1:2        
        if strcmp(stmtype, 'rc')
            prefidx = size(serdata.session(i).results.cond(k).lfpstm.lfp_stm.mean, 1) - 1;
        else
            [~, prefidx] = max(serdata.session(i).results.cond(k).lfpstm.stm.tu{1}.mean);
        end
        para.cond(k).lfpstm.lfp.trace(i, :) = ...
            serdata.session(i).results.cond(k).lfpstm.lfp_stm.mean(1+prefidx, 1:len_trace_lfp);
        for u = 1:3
            para.cond(k).lfpstm.lfp.period(u).freq(i,:) = ...
                serdata.session(i).results.cond(k).lfpstm.period(u).freq(1+prefidx,:);
            para.cond(k).lfpstm.lfp.period(u).power(i,:) = ...
                serdata.session(i).results.cond(k).lfpstm.period(u).pow_avg(1+prefidx, :);
        end
    end
end

% LFP trace and power ==================
plot_lfp_trace(para)
set(gcf, 'Name', [stmtype ': LFP trace & power'], 'NumberTitle','off')

% LFP tuning to stimuli =================
stm = [];
bands = serdata.session(1).results.cond(1).lfpstm.bands;
para.bands = bands;
for i = 1:lenr
    stmidx = serdata.session(i).results.cond(1).lfpstm.stm.vals < 1000;
    stm = [stm, serdata.session(i).results.cond(1).lfpstm.stm.vals(stmidx)];
    for k = 1:2       
        for b = 1:5
            for u = 1:3
                try
                    para.cond(k).lfpstm.power.(bands{b})(u,i) = ...
                        nanmean(serdata.session(i).results.cond(k).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx));
                catch
                    para.cond(k).lfpstm.power.(bands{b})(u,i) = nan;
                end
            end
        end
    end
end
[stm, counts] = uniquecount(stm);
stm(counts < mean(counts)*0.3) = [];
lenstm = length(stm);
para.stm = stm;
for k = 1:2
    for u = 1:3
        para.cond(k).lfpstm.period(u).power_tuning.delta = nan(lenr, lenstm);
        para.cond(k).lfpstm.period(u).power_tuning.theta = nan(lenr, lenstm);
        para.cond(k).lfpstm.period(u).power_tuning.alpha = nan(lenr, lenstm);
        para.cond(k).lfpstm.period(u).power_tuning.beta = nan(lenr, lenstm);
        para.cond(k).lfpstm.period(u).power_tuning.gamma = nan(lenr, lenstm);
    end
end
for i = 1:lenr
   stmidx = serdata.session(i).results.cond(1).lfpstm.stm.vals < 1000;
   idx = ismember(serdata.session(i).results.cond(1).lfpstm.stm.vals(stmidx), stm);
   idx2 = ismember(stm, serdata.session(i).results.cond(1).lfpstm.stm.vals(stmidx));
    for b = 1:5
        for u = 1:3
            try
                pow = serdata.session(i).results.cond(1).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx);
                para.cond(1).lfpstm.period(u).power_tuning.(bands{b})(i,idx2) = pow(idx);
                pow = serdata.session(i).results.cond(2).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx);
                para.cond(2).lfpstm.period(u).power_tuning.(bands{b})(i,idx2) = pow(idx);
            catch
                continue
            end

            % normalization
            m = nanmean([para.cond(1).lfpstm.period(u).power_tuning.(bands{b})(i,:), ...
                para.cond(2).lfpstm.period(u).power_tuning.(bands{b})(i,:)]);
            para.cond(1).lfpstm.period(u).power_tuning.(bands{b})(i,:) = ...
                para.cond(1).lfpstm.period(u).power_tuning.(bands{b})(i,:)/m;
            para.cond(2).lfpstm.period(u).power_tuning.(bands{b})(i,:) = ...
                para.cond(2).lfpstm.period(u).power_tuning.(bands{b})(i,:)/m;
        end
    end
end

[h1, h2] = plot_pow_scatter(para);
set(h1, 'Name', [stmtype ': Band power'], 'NumberTitle','off')
if ~strcmp(para.stmtype, 'rc')
    set(h2, 'Name', [stmtype ': Tuning of LFP power'], 'NumberTitle','off')
end

%% subfunctions
function unity_plot(x,y,bi)
if nargin < 3
    bi = ones(length(x), 1);
end
oks = ~isnan(x) & ~isnan(y);
x = x(oks);
y = y(oks);
bi = bi(oks);
minima = min([x, y]); 
maxima = max([x, y]); 
dist = maxima - minima;
minima = minima - dist*0.1;
maxima = maxima + dist*0.1;
plot([minima maxima], [minima maxima], '-', 'color', 0.6*ones(1,3))
hold on;
for i = 1:length(bi)
    if bi(i)==1
        col = [1 0 0];
    else
        col = [0 0 0];
    end
    scatter(x(i), y(i), 20, 'o',...
            'markerfacecolor', col, 'markeredgecolor', col, ...
            'markerfacealpha', 0.4, 'markeredgealpha', 0.8);
        hold on;
end
if nargin < 3 || length(unique(bi))==1
    try
        p2 = signrank(x, y);
    catch
        p2 = nan;
    end
    text(minima + 0.7*dist, minima + 0.15*dist, ['p = ' num2str(p2)])
else
    try
        p1 = signrank(x(bi==0),y(bi==0));
        p2 = signrank(x(bi==1),y(bi==1));
    catch
        p1 = nan;
        p2 = nan;
    end
    text(minima + 0.7*dist, minima + 0.1*dist, ['nacl: p = ' num2str(p1)])
    text(minima + 0.7*dist, minima + 0.3*dist, ['5ht: p = ' num2str(p2)])
end

axis([minima maxima minima maxima])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square

function plot_lfp_trace(para)
figure;
cols = [[0 0 0]; [1 0 0]];
switch para.pairtype
    case {'drug', 'ps_drug'}
        titlenames = {'NaCl', '5HT'};
        for k = 1:2
            for d = 1:2
                % LFP traces
                subplot(2,6,[1:3] + 3*(k - 1))
                me = nanmean(para.cond(d).lfpstm.lfp.trace(para.is5ht==k-1,:), 1);
                x = [1:length(me)]/1000;
                sem = nanstd(para.cond(d).lfpstm.lfp.trace(para.is5ht==k-1,:), [], 1)...
                    /sqrt(sum(para.is5ht==k-1));
                fill_between(x, me - sem, me + sem, cols(d,:))
                hold on;
                plot(x, me, '-', 'color', cols(d,:))
                if d==2
                    % format
                    yy = get(gca, 'YLim');
                    plot([0 0],yy, '-k')
                    xlim([0 (length(me)+1)/1000])
                    ylim(yy)
                    xlabel('time (s)')
                    ylabel('LFP')
                    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
                    
                    % stats
                    pval = nan(1, length(x));
                    for t = 1:length(x)
                        pval(t) = signrank(para.cond(1).lfpstm.lfp.trace(:, t), ...
                            para.cond(2).lfpstm.lfp.trace(:, t));
                    end
                    p = nan(1, t);
                    p(pval < 0.05/t) = 1;
                    plot(x, yy(2)*p, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
                    title(titlenames{k})
                end

                for u = 1:3
                    % power
                    subplot(2,6, u + 3*(k - 1) + 6)                
                    x = nanmean(para.cond(d).lfpstm.lfp.period(u).freq, 1);
                    me = nanmean(para.cond(d).lfpstm.lfp.period(u).power(para.is5ht==k-1,:), 1);
                    sem = nanstd(para.cond(d).lfpstm.lfp.period(u).power(para.is5ht==k-1,:), [], 1)...
                        /sqrt(sum(para.is5ht==k-1));
                    fill_between(x, me - sem, me + sem, cols(d,:))
                    hold on;
                    plot(x, me, '-', 'color', cols(d,:))
                    hold on;
                    if d==2
                        % format
                        title(['period ' num2str(u)])
                        xlabel('frequency (Hz)')
                        ylabel('power')
                        set(gca, 'YScale', 'log')
                        yy = get(gca, 'YLim');
                        ylim(yy)
                        hold on;                        
                        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
                        
                        % stats
                        pval = nan(1, length(x));
                        for t = 1:length(x)
                            pval(t) = signrank(log(para.cond(1).lfpstm.lfp.period(u).power(para.is5ht==k-1,t)),...
                                log(para.cond(2).lfpstm.lfp.period(u).power(para.is5ht==k-1,t)));
                        end
                        p = nan(1, t);
                        p(pval < 0.05/t) = 1;
                        plot(x, yy(2)*p, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
                    end
                end
            end     
        end
        
    otherwise
        titlenames = 'baseline';
        for d = 1:2
            % LFP trace
            subplot(2,3,1:3)
            me = nanmean(para.cond(d).lfpstm.lfp.trace, 1);
            x = [1:length(me)]/1000;
            sem = nanstd(para.cond(d).lfpstm.lfp.trace, [], 1)...
                /sqrt(size(para.cond(d).lfpstm.lfp.trace, 1));
            fill_between(x, me - sem, me + sem, cols(d,:))
            hold on;
            plot(x, me, '-', 'color', cols(d,:))
            hold on;
            if d==2
                % format
                yy = get(gca, 'YLim');
                plot([0 0],yy, '-k')
                xlim([0 (length(me)+1)/1000])
                ylim(yy)
                xlabel('time (s)')
                ylabel('LFP')
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
                title(titlenames)
                
                % stats
                pval = nan(1, length(me));
                for t = 1:length(me)
                    pval(t) = signrank(para.cond(1).lfpstm.lfp.trace(:,t),...
                        para.cond(2).lfpstm.lfp.trace(:,t));
                end
                p = nan(1, t);
                p(pval < 0.05/t) = 1;
                plot(x, yy(2)*p, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
            end

            for u = 1:3
                % power
                subplot(2,3,3+u)                
                x = nanmean(para.cond(d).lfpstm.lfp.period(u).freq, 1);
                me = nanmean(para.cond(d).lfpstm.lfp.period(u).power, 1);
                sem = nanstd(para.cond(d).lfpstm.lfp.period(u).power, [], 1)...
                    /sqrt(size(para.cond(d).lfpstm.lfp.period(u).power, 1));
                fill_between(x, me - sem, me + sem, cols(d,:))
                hold on;
                plot(x, me, '-', 'color', cols(d,:))
                hold on;
                if d==2
                    % format
                    title(['period ' num2str(u)])
                    xlabel('frequency (Hz)')
                    ylabel('power')
                    set(gca, 'YScale', 'log')
                    yy = get(gca, 'YLim');        
                    ylim(yy)
                    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
                    
                    % stats                    
                    pval = nan(1, length(x));
                    for t = 1:length(x)
                        pval(t) = signrank(log(para.cond(1).lfpstm.lfp.period(u).power(:,t)),...
                            log(para.cond(2).lfpstm.lfp.period(u).power(:,t)));
                    end
                    p = nan(1, t);
                    p(pval < 0.05/t) = 1;
                    plot(x, yy(2)*p, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
                end
            end
        end
end

function [h1, h2] = plot_pow_scatter(para)
h1 = nan; h2 = nan;
cols = [[0 0 0]; [1 0 0]];
linestyles = {'-', '--'};
switch para.pairtype
    case {'sc', 'ps'}
        para.is5ht = zeros(1, length(para.is5ht));
        psize = 3;
    otherwise
        psize = 6;
end
for b = 1:5
    for u = 1:3
        h1 = figure(2);
        % effect size of 5HT on tuning
        subplot(3,5,(u-1)*5+b)
        unity_plot(para.cond(1).lfpstm.power.(para.bands{b})(u,:),...
            para.cond(2).lfpstm.power.(para.bands{b})(u,:), para.is5ht)
        if u==1
            title(para.bands{b})
        end
        if b==1
            xlabel(para.names{1})
            ylabel({['period ' num2str(u)], para.names{2}})
        end

        % tuning curve
        if ~strcmp(para.stmtype, 'rc')
            h2 = figure(3);
            for l = 1:psize/3
                for d = 1:2
                    subplot(5, psize, psize*(b-1)+3*(l-1)+u)
                    me = nanmean(para.cond(d).lfpstm.period(u).power_tuning.(para.bands{b})(para.is5ht==l-1,:), 1);
                    sem = nanstd(para.cond(d).lfpstm.period(u).power_tuning.(para.bands{b})(para.is5ht==l-1,:), [], 1)/...
                        sqrt(size(para.cond(d).lfpstm.period(u).power_tuning.(para.bands{b})(para.is5ht==l-1,:),1));
                    errorbar(para.stm, me, sem, 'color', cols(l,:), 'linestyle', linestyles{d}, 'CapSize', 0)
                    hold on;
                end
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
            end            
            if b==1 
                xlabel(para.stmtype)
                ylabel('normalized power')
            end
        end
    end
end

function [stm, counts] = uniquecount(s)
stm = unique(s);
nuni = length(stm);
counts = zeros(1, nuni);
for i = 1:nuni
    counts(i) = sum(s==stm(i));
end