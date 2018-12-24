function visualizer(para, name0, name2, prefix, save_flag)
%%
% visualize results from LFPanalyzer.m

name = {name0, name2};
periodname = {'before stm.', 'after stm', 'sustained'};
if mean(ismember('gpfs0', cd))==1
    savedir = '/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Figures_kk/';
else
    savedir = 'Z:\Katsuhisa\serotonin_project\LFP_project\Figures_kk\';
end

% % spike-triggered LFP =========================
% h = figure;
% wnd = 0.1; % was 0.3
% cutf = 50;
% for k = 1:2
%     switch k
%         case 1
%             col = zeros(1,3);
%         case 2
%             col = [1,0,0];
%     end
%     
%     % stlfp
%     subplot(1,3,1)
%     fill_between(-wnd:0.001:wnd, para.cond(k).lfpstm.period(2).stlfp.avg_stlfp(end,:) - ...
%         para.cond(k).lfpstm.period(2).stlfp.sem_stlfp(end,:), ...
%         para.cond(k).lfpstm.period(2).stlfp.avg_stlfp(end,:) + para.cond(k).lfpstm.period(2).stlfp.sem_stlfp(end,:), col);
%     hold on;
%     plot(-wnd:0.001:wnd, para.cond(k).lfpstm.period(2).stlfp.avg_stlfp(end,:), '-','color',col);
%     hold on;
%     
%     % power and frequency (< 30Hz)
%     freq = para.cond(k).lfpstm.period(2).stlfp.FREQ(end,:);
%     subplot(1,3,2)
%     plot(freq(freq<cutf), para.cond(k).lfpstm.period(2).stlfp.POW(end,freq<cutf), '-','color',col);
%     set(gca, 'YScale', 'log')
%     hold on;
%     
%     % power and frequency (>= 30Hz)
%     subplot(1,3,3)
%     plot(freq(freq>=cutf), para.cond(k).lfpstm.period(2).stlfp.POW(end,freq>=cutf), '-','color',col);
%     set(gca, 'YScale', 'log')
%     hold on;
% end
% 
% % cosmetics
% subplot(1,3,1)
% yy = get(gca, 'YLim');
% hold on;
% plot([0 0], yy, '-k')
% xlabel('time rel:spike [s]');
% ylabel('avg LFP +/- SEM (\muV)');
% xlim([-wnd wnd]);
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
% subplot(1,3,2)
% xlabel('frequency [Hz]');
% ylabel('power');
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% xlim([0 30])
% 
% subplot(1,3,3)
% xlabel('frequency [Hz]');
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% xlim([30 100])
% 
% % stlfp0 = para.cond(1).stlfp.mean([-wnd:0.001:wnd]==0);
% % stlfp2 = para.cond(2).stlfp.mean([-wnd:0.001:wnd]==0);
% 
% figname = [name{1}, 'VS', name{2}, '_' 'spLFP_' prefix];
% set(h, 'Name', figname,'NumberTitle','off')
% if save_flag==1
%     savefig(h, strcat(savedir, figname, '.fig'))
% end

% stimulus type =====================================
h = figure;    
for k = 1:2       
    lenv = length(para.cond(k).lfpstm.stm.vals);
    col = jet(lenv);
    if lenv==1
        if k==1
            lc = 'k';
        else
            lc = 'r';
        end
          
        % LFP traces
        subplot(2,3,1:3)
        plot(para.cond(k).lfpstm.ts, para.cond(k).lfpstm.lfp_stm.mean, '-', 'color', lc)
        hold on;
        yy = get(gca, 'YLim');
        if k==2
            xlabel('time (s)')
        end
        message = sprintf([name{k} '\n LFP']);
        xlim([-0.2 max(para.cond(k).lfpstm.ts(~isnan(para.cond(k).lfpstm.ts)))])
        ylabel(message)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        if k==2
            xx = get(gca, 'XLim');
            yy = get(gca, 'YLim');
            text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.2*(yy(2)-yy(1)), name{1}, 'color', 'k')
            text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.1*(yy(2)-yy(1)), name{2}, 'color', 'r')
        end
        ylim(yy)
        plot([0 0], yy, '-k')
        
        % power
        for u = 1:3
            subplot(2,3,3+u)
            plot(para.cond(k).lfpstm.period(u).freq, para.cond(k).lfpstm.period(u).pow_avg, '-', 'color', lc)
            set(gca, 'YScale', 'log')
            hold on;
            yy = get(gca, 'YLim');
            plot([0 0], yy, '-k')
            ylim(yy)
            if k==2
                xlabel('frequency (Hz)')
            end
            if u==1
                ylabel('power')
            end
            title(periodname{u})
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end
    else % 4 stimuli per trial
        for s = 1:lenv
            % LFP traces
            subplot(2,6,[1:3]+3*(k-1))
            plot(para.cond(k).lfpstm.ts, para.cond(k).lfpstm.lfp_stm.mean(s,:), '-', 'color', col(s,:))
            hold on;

            % power
            for u = 1:3
                subplot(2,6,6+u+3*(k-1))
                plot(para.cond(k).lfpstm.period(u).freq, para.cond(k).lfpstm.period(u).pow_avg(s,:), '-', 'color', col(s,:))
                set(gca, 'YScale', 'log')
                xlabel('frequency (Hz)')
                ylabel('power')
                title(['period ' num2str(u)])
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
                hold on;
            end
        end
        % cosmetics
        subplot(2,6,[1:3]+3*(k-1))
        if k==1
            yy = get(gca, 'YLim');
        end
        plot([0 0], yy, '-k')
        if k==2
            xlabel('time (s)')
        end
        message = sprintf([name{k} '\n LFP']);
        xlim([-0.2 max(para.cond(k).lfpstm.ts(~isnan(para.cond(k).lfpstm.ts)))])
        ylim(yy)
        ylabel(message)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
end

figname = [name{1}, 'VS', name{2}, '_' 'LFPandPOW_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

% stimulus tuning by LFP bands ============
h = figure;    
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
       
    for k = 1:2
        switch k
            case 1
                lstyle = '-ok';
            case 2
                lstyle = '-^r';
        end
        for u = 1:3
            % tuning curve
            subplot(3,5,5*(u-1)+b)
            plot(para.cond(k).lfpstm.stm.vals(para.cond(k).lfpstm.stm.vals < 1000), ...
                para.cond(k).lfpstm.period(u).lfp_stm_wave(b).pow(para.cond(k).lfpstm.stm.vals < 1000), ...
                lstyle)   
            hold on;
            if b==1 && k==2
                ylabel({['period ' num2str(u)], 'power'})
                xlabel(para.cond(k).lfpstm.stm.param)
            end
            if u==1 && k==2
                title(band)
            end
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end
    end
end

figname = [name{1}, 'VS', name{2}, '_' 'LFPtuning_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

% encoding ============
h = figure;
indnames = {'spike count', 'LFP (uV)', 'delta', 'theta', 'alpha', 'beta', 'gamma'};
fnames = {'reliability', 'selectivity', 'snr2', 'discriminability', 'metabcost'};
leni = length(indnames);
lenf = length(fnames);
cols = {[0 0 0], [1 0 0]};
for i = 1:leni
    for k = 1:2
        % tuning
        subplot(lenf+1, leni, i)
        stm = para.cond(k).lfpstm.stm.tu{i}.unistm;
        me = para.cond(k).lfpstm.stm.tu{i}.mean;
        sem = para.cond(k).lfpstm.stm.tu{i}.std./sqrt(para.cond(k).lfpstm.stm.tu{i}.ntr);
        errorbar(stm, me, sem, '-', 'color', cols{k}, 'capsize', 0)
        hold on;
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        if k==2
            title(indnames{i})
            if i==1
                ylabel('tuning')
            end
        end
        
        % other indices
        for j = 1:lenf
            subplot(lenf+1, leni, i+j*leni)
            v = para.cond(k).lfpstm.stm.tu{i}.(fnames{j});
            if length(v) > 1
                [~, prefidx] = max(para.cond(k).lfpstm.stm.tu{1}.mean);
                v = v(prefidx);
            end
            bar(k, v, 'k')
            hold on;          
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if j==lenf && k==2
                set(gca, 'XTick', [1 2], 'XTickLabel', {name0, name2})
            else
                set(gca, 'XTick', [1 2], 'XTickLabel', {'', ''})
            end
            if i==1 && k==2
                ylabel(fnames{j})
            end
        end
    end
end
figname = [name{1}, 'VS', name{2}, '_' 'encoding_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

% 
% 
% % stimulus tuning by stLFP bands ============
% h = figure;    
% 
% for b = 1:5
%     switch b
%         case 1
%             band = 'delta';
%         case 2
%             band = 'theta';
%         case 3
%             band = 'alpha';
%         case 4
%             band = 'beta';
%         case 5
%             band = 'gamma';
%     end
%     subplot(1,5,b)
%     yy = get(gca, 'YLim');
%     plot([0 0], yy, '-k')    
%     
%     % tuning curve
%     subplot(1,5,b)
%     plot(para.cond(1).lfpstm.stm.vals(para.cond(1).lfpstm.stm.vals < 1000), ...
%         para.cond(1).lfpstm.period(2).stlfp.band(para.cond(1).lfpstm.stm.vals < 1000, b), ...
%         '-ok')
%     hold on;
%     plot(para.cond(2).lfpstm.stm.vals(para.cond(2).lfpstm.stm.vals < 1000), ...
%         para.cond(2).lfpstm.period(2).stlfp.band(para.cond(2).lfpstm.stm.vals < 1000, b), ...
%         '--^r')
%     if b==1
%         ylabel('power')
%         xlabel(para.cond(1).lfpstm.stm.param)
%         title(band)
%     end
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% end
% 
% figname = [name{1}, 'VS', name{2}, '_' 'stLFPtuning_' prefix];
% set(h, 'Name', figname,'NumberTitle','off')
% if save_flag==1
%     savefig(h, strcat(savedir, figname, '.fig'))
% end


% % Spectrogram ===================== 
% h = figure;
% % crange = [0, 1];
% pos = 1;
% for k = 1:2        
%     for s = 1:lenv        
%         for u = 1:3
%             subplot(3, 3*lenv, pos)
%             plot_matrix(para.cond(k).lfpstm.period(u).spectrogram.S{s}, ...
%                 para.cond(k).lfpstm.period(u).spectrogram.t{s}, ...
%                 para.cond(k).lfpstm.period(u).spectrogram.f{s})
% %             imagesc(para.cond(k).lfpstm.period(u).spectrogram.t{s}, ...
% %                 para.cond(k).lfpstm.period(u).spectrogram.f{s},...
% %                 mag2db(para.cond(k).lfpstm.period(u).spectrogram.S{s})')
%             pos = pos + 1;
%             axis xy;
%             colormap('jet')
% %             c = caxis;
%             colorbar('off')            
% %             if crange(1) > c(1)
% %                 crange(1) = c(1);
% %             end
% %             if crange(2) < c(2)
% %                 crange(2) = c(2);
% %             end
% %             caxis(crange)
%             if k==1 && s==1 && u==1
%                 title({[para.cond(k).lfpstm.stm.param, ' = '], ...
%                     num2str(para.cond(k).lfpstm.stm.vals(s))})
%             elseif u == 1 && k == 1
%                 title(num2str(para.cond(k).lfpstm.stm.vals(s)))
%             else
%                 title('')
%             end        
%             if s==1
%                 message = sprintf([name{k} ' \n frequency (Hz)']);
%                 ylabel(message);
%             else
%                 ylabel('')
%             end
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%         end        
%     end
% end
% 
% % difference
% % crange = [-1 1];
% for s = 1:lenv
%     for u = 1:3
%         subplot(3, 3*lenv, pos)
%         imagesc(para.cond(1).lfpstm.period(u).spectrogram.t{s}, ...
%                 para.cond(2).lfpstm.period(u).spectrogram.f{s},...
%                 (mag2db(para.cond(1).lfpstm.period(u).spectrogram.S{s})...
%                 - mag2db(para.cond(2).lfpstm.period(u).spectrogram.S{s}))')
%         pos = pos + 1;
%         axis xy;
%         colormap('jet')
%         colorbar('off')
% 
% %         c = caxis;
% %         if crange(1) > c(1)
% %             crange(1) = c(1);
% %         end
% %         if crange(2) < c(2)
% %             crange(2) = c(2);
% %         end
% %         caxis(crange)
%         title('')
%         if s==1
%             xlabel('time [s]');
%             message = sprintf([name{1} ' - ' name{2} ' \n frequency (Hz)']);
%             ylabel(message);
%         else
%             xlabel('')
%             ylabel('')
%         end
%         set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%     end
% end
% 
% figname = [name{1}, 'VS', name{2}, '_' 'Spectrogram_' prefix];
% set(h, 'Name', figname,'NumberTitle','off')
% if save_flag==1
%     savefig(h, strcat(savedir, figname, '.fig'))
% end

%  % Spike-LFP coherency ===================== 
% h = figure;
% yrange1 = [0 0];
% yrange2 = [0 0];
% for s = 1:lenv
%     subplot(2, lenv, s)
%     plot(para.cond(1).lfpstm.period(2).coherence.f{s}, para.cond(1).lfpstm.period(2).coherence.C{s}, '-', 'color', col(s,:))
%     hold on;
%     plot(para.cond(2).lfpstm.period(2).coherence.f{s}, para.cond(2).lfpstm.period(2).coherence.C{s}, '--', 'color', col(s,:))
%     title([para.cond(1).lfpstm.stm.param, ' = ' ...
%             num2str(para.cond(1).lfpstm.stm.vals(s))])
%     yy = get(gca, 'YLim');
%     if yrange1(1) > yy(1)
%         yrange1(1) = yy(1);
%     end
%     if yrange1(2) < yy(2)
%         yrange1(2) = yy(2);
%     end
%     if s==1
%         ylabel('coherence');
%     end
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%     
%     subplot(2, lenv, s + lenv)
%     plot(para.cond(1).lfpstm.period(2).coherence.f{s}, ...
%         zeros(1, length(para.cond(1).lfpstm.period(2).coherence.f{s})), '-k')
%     hold on;
%     plot(para.cond(1).lfpstm.period(2).coherence.f{s}, para.cond(1).lfpstm.period(2).coherence.C{s} - ...
%         para.cond(2).lfpstm.period(2).coherence.C{s}, ':', 'color', col(s,:))
%     yy = get(gca, 'YLim');
%     if yrange2(1) > yy(1)
%         yrange2(1) = yy(1);
%     end
%     if yrange2(2) < yy(2)
%         yrange2(2) = yy(2);
%     end
%     if s==1
%         xlabel('frequency (Hz)');
%         ylabel('\Delta coherence')
%     end
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% end
% for s = 1:lenv
%     subplot(2, lenv, s)
%     ylim(yrange1)
%     
%     subplot(2, lenv, s + lenv)
%     ylim(yrange2)
% end
% 
% figname = [name{1}, 'VS', name{2}, '_' 'SpkLfpCoherence_' prefix];
% set(h, 'Name', figname,'NumberTitle','off')
% if save_flag==1
%     savefig(h, strcat(savedir, figname, '.fig'))
% end

%  % Spike-LFP phase ===================== 
% h = figure;
% for i = 1:lenv
%     try
%         subplot(2, lenv, i)
%         polarhistogram(para.cond(1).lfpstm.period(2).coherence.phi{i}, ...
%             length(para.cond(1).lfpstm.period(2).coherence.phi{i}), ...
%             'FaceColor','red','FaceAlpha',.3);
%         title([para.cond(1).lfpstm.stm.param, ' = ' ...
%                num2str(para.cond(1).lfpstm.stm.vals(i))])
%     catch
%         disp(['stimulus index ' num2str(i) ' error'])
%     end
%     try
%         subplot(2, lenv, i+lenv)
%         polarhistogram(para.cond(2).lfpstm.period(2).coherence.phi{i}, ...
%             length(para.cond(2).lfpstm.period(2).coherence.phi{i}), ...
%             'FaceColor','red','FaceAlpha',.3);
%     catch
%         disp(['stimulus index ' num2str(i) ' error'])
%     end
% end
% 
% figname = [name{1}, 'VS', name{2}, '_' 'LfpPhaseOfSpikes_' prefix];
% set(h, 'Name', figname,'NumberTitle','off')
% if save_flag==1
%     savefig(h, strcat(savedir, figname, '.fig'))
% end