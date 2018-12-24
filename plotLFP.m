function plotLFP(LFPinfo, stmtype)
%% plot LFP data across sessions
%
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\LFPinfo.mat')
%
% written by Katsuhisa (19.03.18)
% this is old
% ++++++++++++++++++++++++++++++++

close all;
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_useufulfunc'))

%% LFP data: baseline vs drug
% stimulus type
row = [];
is5ht = [];
for i = 1:length(LFPinfo.session)    
    if LFPinfo.session(i).exist==1 && LFPinfo.session(i).goodunit==1 ...
        && isstruct(LFPinfo.session(i).results) && strcmp(stmtype, LFPinfo.session(i).results.stimulus)  
        row = [row, i];        
        if strcmp(LFPinfo.session(i).results.drugname, '5HT')
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
    

% wnd = 0.1;  % was 0.3
for f = 1:3
    switch f
        case 1
            fieldname = 'drug';
            xname = 'baseline';
            yname = 'drug';
        case 2
            fieldname = 'ps_base';
            xname = 'S-ps';
            yname = 'L-ps';
        case 3
            fieldname = 'ps_drug';
            xname = 'S-ps';
            yname = 'L-ps';
    end
    
    % data extraction
    if strcmp(stmtype, 'rc')
        len_trace_lfp = 2000;
    else
        len_trace_lfp = 450;
    end
%     len_trace_lfp = size(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.lfp_stm.mean, 2);
%     len_trace_lfp = 2000;
    xfreq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.f;
    len_pow_lfp = length(xfreq);
    for k = 1:2
%         para.cond(k).lfpstm.stlfp.amp = nan(1, lenr);
%         para.cond(k).lfpstm.stlfp.trace = nan(lenr, len_trace_sta);
%         para.cond(k).lfpstm.stlfp.power = nan(lenr, len_pow_sta);
        para.cond(k).lfpstm.lfp.trace = nan(lenr, len_trace_lfp);
        para.cond(k).lfpstm.lfp.power = nan(lenr, len_pow_lfp);
    end
    for i = 1:lenr
        for k = 1:2            
%             para.cond(k).lfpstm.stlfp.trace(i,:) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.avg_stlfp(end, :);
%             para.cond(k).lfpstm.stlfp.pow(i,:) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.pow(end,:)';
%             para.cond(k).lfpstm.stlfp.amp(i) = ...
%                 min(para.cond(k).lfpstm.stlfp.trace(i,:));
            para.cond(k).lfpstm.lfp.trace(i, :) = ...
                LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.lfp_stm.mean(end,1:len_trace_lfp);
            for u = 1:3
                para.cond(k).lfpstm.lfp.period(u).freq(i,:) = ...
                    LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.period(u).freq(end,:);
                para.cond(k).lfpstm.lfp.period(u).power(i,:) = ...
                    LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.period(u).pow_avg(end,:);
            end
        end
    end
    
%     % stLFP ==========================
%     h = figure;
%     len_trace_sta = size(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stlfp.avg_stlfp, 2);
%     len_pow_sta = size(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stlfp.pow, 2);
% %     len_trace_lfp = size(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.lfp_stm.mean, 2);
%     len_trace_lfp = 2000;
%     xfreq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.f;
%     len_pow_lfp = length(xfreq);
%     for k = 1:2
%         para.cond(k).lfpstm.stlfp.amp = nan(1, lenr);
%         para.cond(k).lfpstm.stlfp.trace = nan(lenr, len_trace_sta);
%         para.cond(k).lfpstm.stlfp.power = nan(lenr, len_pow_sta);
%         para.cond(k).lfpstm.lfp.trace = nan(lenr, len_trace_lfp);
%         para.cond(k).lfpstm.lfp.power = nan(lenr, len_pow_lfp);
%     end
%     for i = 1:lenr
%         for k = 1:2            
%             para.cond(k).lfpstm.stlfp.trace(i,:) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.avg_stlfp(end, :);
%             para.cond(k).lfpstm.stlfp.pow(i,:) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.pow(end,:)';
%             para.cond(k).lfpstm.stlfp.amp(i) = ...
%                 min(para.cond(k).lfpstm.stlfp.trace(i,:));
%             para.cond(k).lfpstm.lfp.trace(i, :) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.lfp_stm.mean(end,1:2000);
%             para.cond(k).lfpstm.lfp.power(i,:) = ...
%                 LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.pow_stm.mean(end,:);
%         end
%     end
%     if f==1
%         for k = 1:2
%             switch k
%                 case 1
%                     drugname = 'NaCl';
%                 case 2
%                     drugname = '5HT';
%             end
% 
%             % stLFP traces --- 5HT
%             subplot(3,3,1+3*(k-1))
%             me = nanmean(para.cond(1).lfpstm.stlfp.trace(is5ht==k-1,:), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.trace(is5ht==k-1,:), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(-wnd:0.001:wnd, me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(-wnd:0.001:wnd, me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.trace(is5ht==k-1,:), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.trace(is5ht==k-1,:), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(-wnd:0.001:wnd, me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(-wnd:0.001:wnd, me, '-r')
%             hold on;
%             yy = get(gca, 'YLim');
%             plot([0 0],yy, '-k')
%             xlim([-wnd wnd])
%             ylim(yy)
%             xlabel('time (s)')
%             ylabel('LFP')
%             title(drugname)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
% 
%             % stLFP power (<30Hz)
%             subplot(3,3,2+3*(k-1))
%             freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stlfp.freq(end,:);
%             me = nanmean(para.cond(1).lfpstm.stlfp.pow(is5ht==k-1,freq<30), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.pow(is5ht==k-1,freq<30), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(freq(freq<30), me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(freq(freq<30), me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.pow(is5ht==k-1,freq<30), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.pow(is5ht==k-1,freq<30), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(freq(freq<30), me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq(freq<30), me, '-r')
%             hold on;
%             xlim([min(freq(freq<30)) max(freq(freq<30))])
%             xlabel('frequency (Hz)')
%             ylabel('power')
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
% 
%             % stLFP power (>30Hz)
%             subplot(3,3,3+3*(k-1))
%              me = nanmean(para.cond(1).lfpstm.stlfp.pow(is5ht==k-1,freq>30), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.pow(is5ht==k-1,freq>30), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(freq(freq>30), me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(freq(freq>30), me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.pow(is5ht==k-1,freq>30), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.pow(is5ht==k-1,freq>30), [], 1)...
%                 /sqrt(sum(is5ht==k-1));
%             fill_between(freq(freq>30), me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq(freq>30), me, '-r')
%             hold on;
%             xlim([min(freq(freq>30)) max(freq(freq>30))])
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
% 
%             % stLFP scatter
%             subplot(3,3,7)
%             unity_plot(para.cond(1).lfpstm.stlfp.amp, para.cond(2).lfpstm.stlfp.amp, is5ht)
%             title('stLFP')
%             xlabel(xname)
%             ylabel(yname)
%         end
%         set(h, 'Name', [fieldname ': stLFP'], 'NumberTitle','off')
%         
%     else % pupil size modulation
%             % stLFP traces --- 5HT
%             subplot(1,4,1)
%             me = nanmean(para.cond(1).lfpstm.stlfp.trace, 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.trace, [], 1)...
%                 /sqrt(lenr);
%             fill_between(-wnd:0.001:wnd, me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(-wnd:0.001:wnd, me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.trace, 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.trace, [], 1)...
%                 /sqrt(lenr);
%             fill_between(-wnd:0.001:wnd, me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(-wnd:0.001:wnd, me, '-r')
%             hold on;
%             yy = get(gca, 'YLim');
%             plot([0 0],yy, '-k')
%             xlim([-wnd wnd])
%             ylim(yy)
%             xlabel('time (s)')
%             ylabel('LFP')
%             title(drugname)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
% 
%             % stLFP power (<30Hz)
%             subplot(1,4,2)
%             freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stlfp.freq(end,:);
%             me = nanmean(para.cond(1).lfpstm.stlfp.pow(:,freq<30), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.pow(:,freq<30), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq(freq<30), me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(freq(freq<30), me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.pow(:,freq<30), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.pow(:,freq<30), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq(freq<30), me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq(freq<30), me, '-r')
%             hold on;
%             xlim([min(freq(freq<30)) max(freq(freq<30))])
%             xlabel('frequency (Hz)')
%             ylabel('power')
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
% 
%             % stLFP power (>30Hz)
%             subplot(1,4,3)
%              me = nanmean(para.cond(1).lfpstm.stlfp.pow(:,freq>30), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfp.pow(:,freq>30), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq(freq>30), me - sem, me + sem, zeros(1,3))
%             hold on;
%             plot(freq(freq>30), me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfp.pow(:,freq>30), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfp.pow(:,freq>30), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq(freq>30), me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq(freq>30), me, '-r')
%             hold on;
%             xlim([min(freq(freq>30)) max(freq(freq>30))])
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
% 
%             % stLFP scatter
%             subplot(1,4,4)
%             unity_plot(para.cond(1).lfpstm.stlfp.amp, para.cond(2).lfpstm.stlfp.amp)
%             title('stLFP')
%             xlabel(xname)
%             ylabel(yname)
%             set(h, 'Name', [fieldname ': stLFP'], 'NumberTitle','off')
%      end

    % LFP trace and power
    h = figure;
    if f==1
        for k = 1:2
            switch k
                case 1
                    drugname = 'NaCl';
                case 2
                    drugname = '5HT';
            end
            
            subplot(2,4,1+(k-1)*4)
            me = nanmean(para.cond(1).lfpstm.lfp.trace(is5ht==k-1,:), 1);
            sem = nanstd(para.cond(1).lfpstm.lfp.trace(is5ht==k-1,:), [], 1)...
                /sqrt(sum(is5ht==k-1));
            fill_between([1:length(me)]/1000, me - sem, me + sem, zeros(1,3))
            hold on;
            plot([1:length(me)]/1000, me, '-k')
            hold on;
            me = nanmean(para.cond(2).lfpstm.lfp.trace(is5ht==k-1,:), 1);
            sem = nanstd(para.cond(2).lfpstm.lfp.trace(is5ht==k-1,:), [], 1)...
                /sqrt(sum(is5ht==k-1));
            fill_between([1:length(me)]/1000, me - sem, me + sem, [1 0 0])
            hold on;
            plot([1:length(me)]/1000, me, '-r')
            hold on;
            yy = get(gca, 'YLim');
            plot([0 0],yy, '-k')
            xlim([0 (length(me)+1)/1000])
            ylim(yy)
            xlabel('time (s)')
            ylabel('LFP')
            title(drugname)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
            
            for u = 1:3
                subplot(2,4,1+(k-1)*4+u)                
                for c = 1:2
                    switch c
                        case 1
                            col = [0 0 0];
                        case 2
                            col = [1 0 0];
                    end
                    x = nanmean(para.cond(k).lfpstm.lfp.period(u).freq, 1);
                    me = nanmean(para.cond(c).lfpstm.lfp.period(u).power(is5ht==k-1,:), 1);
                    sem = nanstd(para.cond(c).lfpstm.lfp.period(u).power(is5ht==k-1,:), [], 1)...
                        /sqrt(sum(is5ht==k-1));
                    fill_between(x, me - sem, me + sem, col)
                    hold on;
                    plot(x, me, '-', 'color', col)
                    hold on;
                end
                p = nan(1, length(x));
                for t = 1:length(x)
                    p(t) = signrank(para.cond(1).lfpstm.lfp.period(u).power(is5ht==k-1,t),...
                        para.cond(2).lfpstm.lfp.period(u).power(is5ht==k-1,t));
                end
                pval = nan(1, t);
                pval(p < 0.05/t) = 1;
                title(['period ' num2str(u)])
                xlabel('frequency (Hz)')
                ylabel('power')
                set(gca, 'YScale', 'log')
                yy = get(gca, 'YLim');
                hold on;
                plot(x, yy(2)*pval, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
            end            
        end
        
    else
            subplot(1,4,1)
            me = nanmean(para.cond(1).lfpstm.lfp.trace, 1);
            sem = nanstd(para.cond(1).lfpstm.lfp.trace, [], 1)...
                /sqrt(lenr);
            fill_between([1:length(me)]/1000, me - sem, me + sem, zeros(1,3))
            hold on;
            plot([1:length(me)]/1000, me, '-k')
            hold on;
            me = nanmean(para.cond(2).lfpstm.lfp.trace, 1);
            sem = nanstd(para.cond(2).lfpstm.lfp.trace, [], 1)...
                /sqrt(sum(is5ht==k-1));
            fill_between([1:length(me)]/1000, me - sem, me + sem, [1 0 0])
            hold on;
            plot([1:length(me)]/1000, me, '-r')
            hold on;
            yy = get(gca, 'YLim');
            plot([0 0],yy, '-k')
            xlim([0 length(me)+1])
            ylim(yy)
            xlabel('time (s)')
            ylabel('LFP')
            title(drugname)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
            
            for u = 1:3
                subplot(1,4,1+u)
                for c = 1:2
                    switch c
                        case 1
                            col = [0 0 0];
                        case 2
                            col = [1 0 0];
                    end
                    x = mean(para.cond(k).lfpstm.lfp.period(u).freq, 1);
                    me = nanmean(para.cond(c).lfpstm.lfp.period(u).power(is5ht==k-1,:), 1);
                    sem = nanstd(para.cond(c).lfpstm.lfp.period(u).power(is5ht==k-1,:), [], 1)...
                        /sqrt(sum(is5ht==k-1));
                    fill_between(x, me - sem, me + sem, col)
                    hold on;
                    plot(x, me, '-', 'color', col)
                    hold on;
                end
                p = nan(1, length(x));
                for t = 1:length(x)
                    p(t) = signrank(para.cond(1).lfpstm.lfp.period(u).power(is5ht==k-1,t),...
                        para.cond(2).lfpstm.lfp.period(u).power(is5ht==k-1,t));
                end
                pval = nan(1, t);
                pval(p < 0.05/t) = 1;
                title(['period ' num2str(u)])
                xlabel('frequency (Hz)')
                ylabel('power')
                set(gca, 'YScale', 'log')
                yy = get(gca, 'YLim');
                hold on;
                plot(x, 1.05*yy(2)*pval, '-', 'color', 0.6*[1 1 1], 'linewidth', 1.5)
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
            end            
            
    end
    set(h, 'Name', [fieldname ': LFP driven by stimulus'], 'NumberTitle','off')

    % LFP tuning to stimuli =================
    stm = [];
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];
        for k = 1:2       
            for b = 1:5
                switch b
                    case 1
                        band = 'delta';
                    case 2
                        band = 'theta';
                    case 3
                        band = 'alpha';
                    case 4
                        band = 'beta';
                    case 5
                        band = 'gamma';
                end
                for u = 1:3
                    try
                        para.cond(k).lfpstm.power.(band)(u,i) = ...
                            nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx));
                    catch
                        para.cond(k).lfpstm.power.(band)(u,i) = nan;
                    end
                end
            end
        end
    end
    [stm, counts] = uniquecount(stm);
    stm(counts < mean(counts)*0.3) = [];
    lenstm = length(stm);
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
       stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
       idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
       idx2 = ismember(stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx));
        for b = 1:5
            switch b
                case 1
                    band = 'delta';
                case 2
                    band = 'theta';
                case 3
                    band = 'alpha';
                case 4
                    band = 'beta';
                case 5
                    band = 'gamma';
            end

            for u = 1:3
                try
                    pow = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx);
                    para.cond(1).lfpstm.period(u).power_tuning.(band)(i,idx2) = pow(idx);
                    pow = LFPinfo.session(row(i)).results.(fieldname).cond(2).lfpstm.period(u).lfp_stm_wave(b).pow(stmidx);
                    para.cond(2).lfpstm.period(u).power_tuning.(band)(i,idx2) = pow(idx);
                catch
                    continue
                end

                % normalization
                m = nanmean([para.cond(1).lfpstm.period(u).power_tuning.(band)(i,:), ...
                    para.cond(2).lfpstm.period(u).power_tuning.(band)(i,:)]);
                para.cond(1).lfpstm.period(u).power_tuning.(band)(i,:) = ...
                    para.cond(1).lfpstm.period(u).power_tuning.(band)(i,:)/m;
                para.cond(2).lfpstm.period(u).power_tuning.(band)(i,:) = ...
                    para.cond(2).lfpstm.period(u).power_tuning.(band)(i,:)/m;
            end
        end
    end
    for b = 1:5
        switch b
            case 1
                band = 'delta';
            case 2
                band = 'theta';
            case 3
                band = 'alpha';
            case 4
                band = 'beta';
            case 5
                band = 'gamma';
        end
        for u = 1:3
            % effect size of 5HT on tuning
            h1 = figure(12+f);
            subplot(3,5,(u-1)*5+b)
            if f==1
                unity_plot(para.cond(1).lfpstm.power.(band)(u,:),...
                    para.cond(2).lfpstm.power.(band)(u,:), is5ht)
            else
                unity_plot(para.cond(1).lfpstm.power.(band)(u,:),...
                    para.cond(2).lfpstm.power.(band)(u,:))
            end
            if u==1
                title(band)
            end
            if b==1
                xlabel(xname)
                ylabel({['period ' num2str(u)], yname})
            end

            % tuning curve
            for l = 1:2
                switch l
                    case 1
                        col = [0 0 0];
                    case 2
                        col = [1 0 0];
                end
                h2 = figure(123+f);
                subplot(5,6,6*(b-1)+3*(l-1)+u)
                me = nanmean(para.cond(1).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:), 1);
                sem = nanstd(para.cond(1).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:), [], 1)/...
                    sqrt(size(para.cond(1).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:),1));
                errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
                hold on;
                me = nanmean(para.cond(2).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:), 1);
                sem = nanstd(para.cond(2).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:), [], 1)/...
                    sqrt(size(para.cond(2).lfpstm.period(u).power_tuning.(band)(is5ht==l-1,:),1));
                errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
            end
        end
        if b==1 
            xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
            ylabel('normalized power')
        end
    end
    set(h1, 'Name', [fieldname ': Band power'], 'NumberTitle','off')
    set(h2, 'Name', [fieldname ': Tuning of LFP power'], 'NumberTitle','off')

%     % stLFP tuning to stimuli =================
%     h = figure;
%     for k = 1:2
%         para.cond(k).lfpstm.stlfppower.delta = nan(1, lenr);
%         para.cond(k).lfpstm.stlfppower.theta = nan(1, lenr);
%         para.cond(k).lfpstm.stlfppower.alpha = nan(1, lenr);
%         para.cond(k).lfpstm.stlfppower.beta = nan(1, lenr);
%         para.cond(k).lfpstm.stlfppower.gamma = nan(1, lenr);
%     end
%     stm = [];
%     for i = 1:lenr
%         stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
%         stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];
% 
%         for k = 1:2    
%             para.cond(k).lfpstm.stlfppower.delta(i) = ...
%                 nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,1), 1);
%             para.cond(k).lfpstm.stlfppower.theta(i) = ...
%                 nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,2), 1);
%             para.cond(k).lfpstm.stlfppower.alpha(i) = ...
%                 nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,3), 1);
%             para.cond(k).lfpstm.stlfppower.beta(i) = ...
%                 nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,4), 1);
%             para.cond(k).lfpstm.stlfppower.gamma(i) = ...
%                 nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,5), 1);
%         end
%     end
%     [stm, counts] = uniquecount(stm);
%     stm(counts < mean(counts)*0.3) = [];
%     lenstm = length(stm);
%     for k = 1:2
%         para.cond(k).lfpstm.stlfppower_tuning.delta = nan(lenr, lenstm);
%         para.cond(k).lfpstm.stlfppower_tuning.theta = nan(lenr, lenstm);
%         para.cond(k).lfpstm.stlfppower_tuning.alpha = nan(lenr, lenstm);
%         para.cond(k).lfpstm.stlfppower_tuning.beta = nan(lenr, lenstm);
%         para.cond(k).lfpstm.stlfppower_tuning.gamma = nan(lenr, lenstm);
%     end
%     for i = 1:lenr
%         stmidx = find(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000);
%         idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
%         idx2 = ismember(stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx));
% 
%         for k = 1:2    
%             bandtu = LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,:);
%             para.cond(k).lfpstm.stlfppower_tuning.delta(i,idx2) = bandtu(idx, 1)';
%             para.cond(k).lfpstm.stlfppower_tuning.theta(i,idx2) = bandtu(idx, 2)';
%             para.cond(k).lfpstm.stlfppower_tuning.alpha(i,idx2) = bandtu(idx, 3)';
%             para.cond(k).lfpstm.stlfppower_tuning.beta(i,idx2) = bandtu(idx, 4)';
%             para.cond(k).lfpstm.stlfppower_tuning.gamma(i,idx2) = bandtu(idx, 5)';
%         end
% 
%         for b = 1:5
%             switch b
%                 case 1
%                     band = 'delta';
%                 case 2
%                     band = 'theta';
%                 case 3
%                     band = 'alpha';
%                 case 4
%                     band = 'beta';
%                 case 5
%                     band = 'gamma';
%            end
% 
%             % normalization
%             m = nanmean([para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:), ...
%                 para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:)]);
%             para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:) = ...
%                 para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:)/m;
%             para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:) = ...
%                 para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:)/m;
%         end
%     end
% 
%     for b = 1:5
%         switch b
%             case 1
%                 band = 'delta';
%             case 2
%                 band = 'theta';
%             case 3
%                 band = 'alpha';
%             case 4
%                 band = 'beta';
%             case 5
%                 band = 'gamma';
%         end
% 
%         % effect size of 5HT on tuning
%         subplot(2,5,b)
%         if f==1
%             unity_plot(para.cond(1).lfpstm.stlfppower.(band),...
%                 para.cond(2).lfpstm.stlfppower.(band), is5ht)
%         else
%             unity_plot(para.cond(1).lfpstm.stlfppower.(band),...
%                 para.cond(2).lfpstm.stlfppower.(band))
%         end
%         title(band)
%         if b==1
%             xlabel(xname)
%             ylabel(yname)
%         end
% 
%         % tuning curve
%         for l = 1:2
%             switch l
%                 case 1
%                     col = [0 0 0];
%                 case 2
%                     col = [1 0 0];
%             end
%             subplot(2,5,b+5)
%             me =nanmean(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), 1);
%             sem = nanstd(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), [], 1)/...
%                 sqrt(size(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:),1));
%             errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), 1);
%             sem = nanstd(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), [], 1)/...
%                 sqrt(size(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:),1));
%             errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
%         end
%         if b==1
%             xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
%             ylabel('normalized power')
%         end
%     end
%     set(h, 'Name', [fieldname ': Tuning of stLFP power'], 'NumberTitle','off')
% 
% 
%     % spike-LFP coherence ==================
%     h = figure;
%     stm = [];
%     for i = 1:lenr
%         stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
%         stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];
%     end
%     [stm, counts] = uniquecount(stm);
%     stm(counts < mean(counts)*0.3) = [];
%     lenstm = length(stm);
%     freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.coherence.f{1};
%     lenf = length(freq);
%     for k = 1:2
%         para.cond(k).lfpstm.coherence_stmval = nan(lenr, lenstm);
%         para.cond(k).lfpstm.coherence_freq = nan(lenr, lenf);
%     end
%     for i = 1:lenr
%         stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
%         idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
%         idxi = find(idx==1);
%         for k = 1:2
%             for s = 1:sum(idx==1)
%                 para.cond(k).lfpstm.coherence_stmval(i, s) = ...
%                     nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.coherence.C{idxi(s)});      
%             end
%             para.cond(k).lfpstm.coherence_freq(i,:) = ...
%                     LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.coherence.C{idxi(s)}';
%         end
%     end
% 
%     if f==1
%         for l = 1:2
%             switch l
%                 case 1
%                     drugname = 'NaCl';
%                 case 2
%                     drugname = '5HT';
%             end
% 
%             % coherence vs frequency
%             subplot(2,2,1+2*(l-1))
%             freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.coherence.f{1};
%             me = nanmean(para.cond(1).lfpstm.coherence_freq(is5ht==l-1,:),1);
%             sem = nanstd(para.cond(1).lfpstm.coherence_freq(is5ht==l-1,:), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq, me - sem, me + sem, [0 0 0])
%             hold on;
%             plot(freq, me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.coherence_freq(is5ht==l-1,:),1);
%             sem = nanstd(para.cond(2).lfpstm.coherence_freq(is5ht==l-1,:), [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq, me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq, me, '-r')
%             hold on;
%             xlabel('frequency (Hz)')
%             ylabel('coherence')
%             title(drugname)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
% 
%             % tuning curve
%             switch l
%                 case 1
%                     col = [0 0 0];
%                 case 2
%                     col = [1 0 0];
%             end
% 
%             subplot(2,2,2+2*(l-1))
%             me = nanmean(para.cond(1).lfpstm.coherence_stmval(is5ht==l-1,:),1);
%             sem = nanstd(para.cond(1).lfpstm.coherence_stmval(is5ht==l-1,:), [], 1)/sqrt(lenr);        
%             errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.coherence_stmval(is5ht==l-1,:),1);
%             sem = nanstd(para.cond(2).lfpstm.coherence_stmval(is5ht==l-1,:), [], 1)/sqrt(lenr);
%             errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
%             xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
%         end
%     else        
%             % coherence vs frequency
%             subplot(1,2,1)
%             freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.coherence.f{1};
%             me = nanmean(para.cond(1).lfpstm.coherence_freq,1);
%             sem = nanstd(para.cond(1).lfpstm.coherence_freq, [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq, me - sem, me + sem, [0 0 0])
%             hold on;
%             plot(freq, me, '-k')
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.coherence_freq,1);
%             sem = nanstd(para.cond(2).lfpstm.coherence_freq, [], 1)...
%                 /sqrt(lenr);
%             fill_between(freq, me - sem, me + sem, [1 0 0])
%             hold on;
%             plot(freq, me, '-r')
%             hold on;
%             xlabel('frequency (Hz)')
%             ylabel('coherence')
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
% 
%             % tuning curve
%             subplot(1,2,2)
%             me = nanmean(para.cond(1).lfpstm.coherence_stmval,1);
%             sem = nanstd(para.cond(1).lfpstm.coherence_stmval, [], 1)/sqrt(lenr);        
%             errorbar(stm, me, sem, 'color', 'k', 'linestyle','-', 'CapSize', 0)
%             hold on;
%             me = nanmean(para.cond(2).lfpstm.coherence_stmval,1);
%             sem = nanstd(para.cond(2).lfpstm.coherence_stmval, [], 1)/sqrt(lenr);
%             errorbar(stm, me, sem, 'color', 'r', 'linestyle','--', 'CapSize', 0)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
%             xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
%      end
% 
% 
%     set(h, 'Name', [fieldname ': spike-LFP coherence'], 'NumberTitle','off')
    
%     % spike-LFP phase ==================
%     h = figure;
%     if f==1
%         phi0 = []; phi1 = []; phi2 = []; phi3 = [];
%         for i = 1:lenr         
%              if is5ht(i)==1
%                 phi0 = [phi0; [LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.coherence.phi{:}]];
%                 phi2 = [phi2; [LFPinfo.session(row(i)).results.(fieldname).cond(2).lfpstm.coherence.phi{:}]];
%              elseif is5ht(i)==0
%                  phi1 = [phi1; [LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.coherence.phi{:}]];
%                  phi3 = [phi3; [LFPinfo.session(row(i)).results.(fieldname).cond(2).lfpstm.coherence.phi{:}]];
%              end
%         end
%         for l = 1:2
%             switch l
%                 case 1
%                     drugname = 'NaCl';
%                     phi_c1 = phi1(:); phi_c2 = phi3(:);
%                 case 2
%                     drugname = '5HT';
%                     phi_c1 = phi0(:); phi_c2 = phi2(:);
%             end
% 
%             % LFP phase of spikes
%             subplot(2,2,1+2*(l-1))
%             polarhistogram(phi_c1, 'FaceColor','red','FaceAlpha',.3);
%             title('baseline')
%             subplot(2,2,2+2*(l-1))
%             polarhistogram(phi_c2, 'FaceColor','red','FaceAlpha',.3);
%             title(drugname)
%         end
%     else
%         phi_c1 = []; phi_c2 = []; 
%         for i = 1:lenr         
%              if is5ht(i)==1
%                 phi_c1 = [phi_c1; [LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.coherence.phi{:}]];
%                 phi_c2 = [phi_c2; [LFPinfo.session(row(i)).results.(fieldname).cond(2).lfpstm.coherence.phi{:}]];
%              end
%         end
%         % LFP phase of spikes
%         subplot(1,2,1)
%         polarhistogram(phi_c1(:), 'FaceColor','red','FaceAlpha',.3);
%         title('S-ps')
%         subplot(1,2,2)
%         polarhistogram(phi_c2(:), 'FaceColor','red','FaceAlpha',.3);
%         title('L-ps')
%     end
% 
%     set(h, 'Name', [fieldname ': spike-LFP phase'], 'NumberTitle','off')
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
if nargin < 3
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

