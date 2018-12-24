function plot_FRonLFP
% check whether firing rate modulated by contrast stimulus affects stLFP
% ampltude
%
% Test whether the magnitude of the spike triggered LFP is depending on the
% spiking activity. We use the data recorded with a 2s stimulus, to avoid
% dominant slow fluctuations.
%
% 
% 04.04.18 Katsuhisa wrote it

close all

%% folder specifications
if mean(ismember('gpfs0', cd))==1
    main_dir = '/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/MagstLFP_vs_contrast/'; 
else
   main_dir = 'Z:\Katsuhisa\serotonin_project\LFP_project\MagstLFP_vs_contrast\';
end

%% load data given by 'check_FRonLFP.m'
if nargin < 1; type = 0; end
load([main_dir 'Data/fr_lfp_sta.mat'])

%% visualization
% stLFP in each session
lens = length([fr_lfp.session]);
wnd = fr_lfp.session(1).results.wnd; % was 0.3
rho = nan(lens, 2);
c = 0;
for i = 1:lens
    if isstruct(fr_lfp.session(i).results)
        c = c + 1;
%         % stLFP amplitude vs co in each unit
%         figure(1);
%         subplot(2,5,c)
%         plot(fr_lfp.session(i).results.stm.vals(2:end), ...
%             fr_lfp.session(i).results.stlfp.peak_stlfp(2:end), ...
%             ':ok')
%         rho(i,1) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
%             fr_lfp.session(i).results.stlfp.peak_stlfp(2:end)', 'type', 'Spearman');
%         if ismember(i, 5)
%             ylabel('peak stLFP amplitude (uV)')
%         end
%         if ismember(i, 10)
%             xlabel('contrast')
%         end
%         title(fr_lfp.session(i).results.ntr)
%         set(gca, 'XScale', 'log')
%         set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
%         % stLFP peak time vs co in each unit
%         figure(2);
%         subplot(2,5,c)
%         plot(fr_lfp.session(i).results.stm.vals(2:end), ...
%             fr_lfp.session(i).results.stlfp.t_peak_stlfp(2:end), ...
%             ':ok')
%         rho(i,2) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
%             fr_lfp.session(i).results.stlfp.t_peak_stlfp(2:end)', 'type', 'Spearman');
%         if ismember(i, 5)
%             ylabel('time of peak stLFP (s)')
%         end
%         if ismember(i, 10)
%             xlabel('contrast')
%         end
%         title(fr_lfp.session(i).results.ntr)
%         set(gca, 'XScale', 'log')
%         set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % fr vs co in each unit
        figure(3);
        subplot(2,5,c)
        errorbar(fr_lfp.session(i).results.stm.tu{1}.unistm, ...
            fr_lfp.session(i).results.stm.tu{1}.mean, ...
            fr_lfp.session(i).results.stm.tu{1}.std./sqrt(fr_lfp.session(i).results.stm.tu{1}.ntr), '-ok')
        rho(i,2) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
            fr_lfp.session(i).results.stm.tu{1}.mean', 'type', 'Spearman');
        if ismember(i, 5)
            ylabel('time of peak stLFP (s)')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
        title([num2str(fr_lfp.session(i).results.ntr) ' trials'])
        set(gca, 'XScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % LFP & stLFP traces
        lenstm = length(fr_lfp.session(i).results.stm.vals(2:end)); 
        col = jet(lenstm);
        pl = nan(1,lenstm);
        leg = cell(1, lenstm);
        for k = 1:lenstm
            % LFP time course
            figure(4);
            subplot(2,5,c)
            pl(k) = plot(fr_lfp.session(i).results.ts, ...
            fr_lfp.session(i).results.lfp_stm.mean(k+1,:), 'color', col(k,:));
            hold on;
            
            % STA 
            figure(5);
            subplot(2,5,c)
            fill_between(-wnd:0.001:wnd, ...
                fr_lfp.session(i).results.stlfp.mean(k+1, :) - fr_lfp.session(i).results.stlfp.sem(k+1, :), ...
                fr_lfp.session(i).results.stlfp.mean(k+1, :) + fr_lfp.session(i).results.stlfp.sem(k+1, :), ...
                col(k,:), 0.4)
            hold on;
            plot(-wnd:0.001:wnd, fr_lfp.session(i).results.stlfp.mean(k+1, :), '-', 'color', col(k,:));

            
            for u = 1:3
                figure(5+u);
                subplot(2,5,c)
                fill_between(-wnd:0.001:wnd, ...
                    fr_lfp.session(i).results.period(u).stlfp.avg_stlfp(k+1, :) - fr_lfp.session(i).results.period(u).stlfp.sem_stlfp(k+1, :), ...
                    fr_lfp.session(i).results.period(u).stlfp.avg_stlfp(k+1, :) + fr_lfp.session(i).results.period(u).stlfp.sem_stlfp(k+1, :), ...
                    col(k,:), 0.4)
                hold on;
                plot(-wnd:0.001:wnd, fr_lfp.session(i).results.period(u).stlfp.avg_stlfp(k+1, :), '-', 'color', col(k,:));
                leg{k} = ['co: ' num2str(fr_lfp.session(i).results.stm.vals(k+1))];
                hold on;
            end
        end
        for ff = 4:8
            figure(ff);
            subplot(2,5,c)
            title([num2str(fr_lfp.session(i).results.ntr) ' trials'])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':k')
            if ff==4
                xx = get(gca, 'XLim');
                for j = 1:length(leg)
                    text(xx(1)+0.5*(xx(2)-xx(1)), yy(1)+0.1*(j-1)*(yy(2)-yy(1)), ...
                        leg{j}, 'color', col(j,:))
                end
%                 legend(pl, leg, 'location', 'southeast')
%                 legend('boxoff')            
                xlim([fr_lfp.session(i).results.ts(1) fr_lfp.session(i).results.ts(end)])
                set(gcf, 'Name', 'LFP traces', 'NumberTitle', 'off')
            else
                xlim([-wnd wnd])
                set(gca, 'XTick', [-wnd 0 wnd])
                set(gcf, 'Name', 'stLFP', 'NumberTitle', 'off')
            end
        end
%         
%         if ismember(i, 5)
%             ylabel('stLFP (uV)')
%         end
%         if ismember(i, 10)
%             xlabel('contrast')
%         end
    end
end
    
% % across units
% figure(4);
% subplot(1,2,1)
% me_hist(rho(:,1))
% xlabel('Spearman rho')
% title('co vs peak stLFP amp.')
% subplot(1,2,2)
% me_hist(rho(:,2))
% xlabel('Spearman rho')
% title('co vs peak stLFP time.')

