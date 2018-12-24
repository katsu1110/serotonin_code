function plot_interaction_lfp(pss, varargin)
% visualize the results from 'interaction_lfp_dataset.m'
% load data here:
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\pss.mat')

if nargin < 1
    load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\pss.mat')
end
if nargin < 2
    celltype_flag = 0;
    celltype = 1;
end

if celltype_flag==1
    % just to get a wave width
    load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\rcdat.mat')
    wdt = arrayfun(@(x) x.wdt(1), dat);
    celltype = 3;
end
close all;

%%
% visualization of results from 'pupil_interaction.m'
pss_orig = pss;
f = 1;
for w = 1:celltype
    if w==1
        prefix = 'all ';
    elseif w==2
        pss.session = pss_orig.session(wdt < 200);
        prefix = 'narrow spk ';
    elseif w==3
        pss.session = pss_orig.session(wdt >= 200);
        prefix = 'broad spk ';
    end
    lenses = length(pss.session);
    for u = 1:2
        switch u
            case 1
                unittype = 'su';
            case 2
                unittype = 'mu';
        end
        % retrieve data
        para.drug(1).cntr.pupiltc = [];
        para.drug(2).cntr.pupiltc = [];
        para.drug(1).drug.pupiltc = [];
        para.drug(2).drug.pupiltc = [];
        para.inter_table = nan(lenses, 4);
        para.perf = nan(lenses, 4);
        para.inter_table_lat = nan(lenses, 4);
        para.type2reg = nan(lenses, 4);
        drugtype = nan(lenses, 1);
        nmin = 10000;
        for s = 1:lenses
            lenps = length(pss.session(s).(['psintr_' unittype]).pupil_timecourse_cntr);
            if lenps < nmin
                nmin = lenps;
            end
        end
        for s = 1:lenses
            if pss.session(s).(['psintr_' unittype '_exist'])==1
                if strcmp(pss.session(s).(['psintr_' unittype]).drugname, '5HT')
                    drugtype(s) = 1;
                    para.drug(2).cntr.pupiltc = [para.drug(2).cntr.pupiltc; ...
                        pss.session(s).(['psintr_' unittype]).pupil_timecourse_cntr(:, 1:nmin)];
                    para.drug(2).drug.pupiltc = [para.drug(2).drug.pupiltc; ...
                        pss.session(s).(['psintr_' unittype]).pupil_timecourse_drug(:, 1:nmin)];
                else
                    drugtype(s) = 0;
                    para.drug(1).cntr.pupiltc = [para.drug(1).cntr.pupiltc; ...
                        pss.session(s).(['psintr_' unittype]).pupil_timecourse_cntr(:, 1:nmin)];
                    para.drug(1).drug.pupiltc = [para.drug(1).drug.pupiltc; ...
                        pss.session(s).(['psintr_' unittype]).pupil_timecourse_drug(:, 1:nmin)];
                end
                para.inter_table(s, :) = reshape(pss.session(s).(['psintr_' unittype]).inter_table, 1, 4);
                para.perf(s, :) = pss.session(s).(['psintr_' unittype]).perf;
                para.inter_table_lat(s,:) = pss.session(s).(['psintr_' unittype]).inter_table_lat;
                para.type2reg(s,[1,3]) = ...
                    [pss.session(s).(['psintr_' unittype]).type2reg.drug(2), pss.session(s).(['psintr_' unittype]).type2reg.base(2)];
                para.type2reg(s,[2,4]) = ...
                    [pss.session(s).(['psintr_' unittype]).type2reg.drug(1), pss.session(s).(['psintr_' unittype]).type2reg.base(1)];
            end
        end
        % visualization        
        for d = 1:2
            figure(f);
            % pupil time-course
            subplot(2,3,1+(d-1)*3)
            me0 = mean(para.drug(d).cntr.pupiltc, 1);
            sem0 = std(para.drug(d).cntr.pupiltc, [], 1)/sqrt(sum(drugtype==d-1));
            me2 = mean(para.drug(d).drug.pupiltc, 1);
            sem2 = std(para.drug(d).drug.pupiltc, [], 1)/sqrt(sum(drugtype==d-1));
            xval = linspace(0, 2000, length(me0));
            fill_between(xval, me0 - sem0, me0 + sem0, [0 0 0])
            hold on;
            plot(xval, me0, '-', 'color', [0 0 0])
            hold on;
            fill_between(xval, me2 - sem2, me2 + sem2, [1 0 0])
            hold on;
            plot(xval, me2, '-', 'color', [1 0 0])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if d==1
                title(['NaCl: N = ' num2str(sum(drugtype==d-1))])
            elseif d==2
                title(['5HT: N = ' num2str(sum(drugtype==d-1))])
            end
            xlim([5 2000])
            ylabel('pupil size (a.u.)')
            xlabel('time (ms)')

            % interaction plot
            subplot(2,3,2+(d-1)*3)
%             me = mean(para.inter_table(drugtype==d-1, :), 1);
%             sem = std(para.inter_table(drugtype==d-1, :), [], 1)/sqrt(sum(drugtype==d-1));
%             errorbar(1:2, me([3,4]), sem([3,4]), '-k', 'capsize', 0)
%             hold on;
%             errorbar(1:2, me([1,2]), sem([1,2]), '-r', 'capsize', 0)
%             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')        
%             xlim([0.7 2.3])
            interaction_plot(para.inter_table(drugtype==d-1, :))
            xlabel({'pupil size', '(small)'})      
            ylabel({'pupil size', '(large)'}) 
%             if d==2
%                 set(gca, 'XTick', 1:2, 'XTickLabel', {'small', 'large'})
%                 xlabel('pupil size')
%             else
%                 set(gca, 'XTick', 1:2, 'XTickLabel', {'', ''})
%             end

            % variance explained
            subplot(2,3,3+(d-1)*3)
            me = mean(para.perf(drugtype==d-1, :), 1);
            sem = std(para.perf(drugtype==d-1, :), [], 1)/sqrt(sum(drugtype==d-1));
            errorbar(1:4, me, sem, '-k', 'capsize', 0)
            xlim([0.7 4.3])
            ylim([0.5 0.75])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if d==2
                set(gca, 'XTick', 1:4, 'XTickLabel', {'stimulus', 'drug', 'pupil', 'drug x pupil'})
                xtickangle(45)
            else
                set(gca, 'XTick', 1:4, 'XTickLabel', {'', '', '', ''})
            end

            figure(f+1);
            % interaction plot (latency)
            subplot(2,4,1+(d-1)*4)
            interaction_plot(para.inter_table_lat(drugtype==d-1, :))
            xlabel({'pupil size', '(small)'})      
            ylabel({'pupil size', '(large)'})
            if d==1
                title('latency')
            end
        
            figure(f+2);
            % interaction plot (type2reg)
            subplot(2,4,1+(d-1)*4)
            interaction_plot(para.type2reg(drugtype==d-1, :), 0)
            xlabel({'gain change'})      
            ylabel({'additive change'})
            if d==1
                title('FR')
            end
            
            set(gcf, 'Name', unittype, 'NumberTitle', 'off')
        end         
        f = f + 3;
    end
end

%     %%
%     % visualization of results from 'pupil_lfp.m'    
%     l = length(pss.session(1).pslfp.interaction);
%     para.ps_lfp_corr = nan(lenses, 2*l);
%     for s = 1:lenses
%         if pss.session(s).pslfp_exist==1
%             para.ps_lfp_corr(s, :) = [pss.session(s).pslfp.corr.control.rho, ...
%                 pss.session(s).pslfp.corr.drug.rho];
%             for c = 1:l
%                 para.interaction(c).table(s, :) = reshape(pss.session(s).pslfp.interaction(c).table, 1, 4);
%             end
%             para.inter_table_lat(s,:) = pss.session(s).pslfp.inter_table_lat;
%             para.type2reg(s, [1,3]) = [pss.session(s).pslfp.type2reg.drug(2), ...
%                 pss.session(s).pslfp.type2reg.base(2)];
%             para.type2reg(s, 3:4) = [pss.session(s).pslfp.type2reg.drug(1), ...
%                 pss.session(s).pslfp.type2reg.base(1)];
%         end
%     end
%     % correlation -----------------------------
%     for d = 1:2        
%         for c = 1:l
%             figure(f);
%             subplot(3,4,c)
% %             if c==12
% %                 a = para.ps_lfp_corr(drugtype==d-1, c+l)
% %                 b = para.ps_lfp_corr(drugtype==d-1, c)
% %             end
%             unity_scatter(para.ps_lfp_corr(drugtype==d-1, c+l), ...
%                 para.ps_lfp_corr(drugtype==d-1, c))
%             set(gcf, 'Name', [prefix ': correlation between ps and lfp'], 'NumberTitle', 'off')
%             switch c
%                 case 1
%                     titlelabel = 'delta';
%                     if d==1
%                          ylabel('NaCl')
%                     elseif d==2
%                         ylabel('5HT')
%                     end
%                 case 2
%                     titlelabel = 'theta';
%                 case 3
%                     titlelabel = 'alpha';
%                 case 4
%                     titlelabel = 'beta';
%                 case 5
%                     titlelabel = 'gamma';
%                 case 6
%                     titlelabel = 'st-delta';
%                 case 7
%                     titlelabel = 'st-theta';
%                 case 8
%                     titlelabel= 'st-alpha';
%                 case 9
%                     titlelabel = 'st-beta';
%                 case 10
%                     titlelabel = 'st-gamma';
%                 case 11
%                     titlelabel = 'LFP';
%                 case 12
%                     titlelabel = 'stLFP amp.';
%             end
%             xlabel('baseline')              
%             title(titlelabel)
%             
%             % interaction ----------------------------------
%             figure(f+1);
%             subplot(3,4,c)
%             v = para.interaction(c).table(drugtype==d-1, :);
%             ok = ~isnan(v(:,1)) | v(:,1)>0;
% %             me = mean(v(ok, :), 1);
% %             sem = std(v(ok, :), [], 1)/sqrt(sum(drugtype==d-1));
% %             errorbar(1:2, me([3,4]), sem([3,4]), '-k', 'capsize', 0)
% %             hold on;
% %             errorbar(1:2, me([1,2]), sem([1,2]), '-r', 'capsize', 0)
% %             set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')        
% %             xlim([0.7 2.3])
%     %         yy = get(gca, 'YLim');
%             interaction_plot(v(ok,:))
%             title(titlelabel)
%             xlabel({'pupil size', '(small)'})      
%             ylabel({'pupil size', '(large)'})   
%         end  
%         set(gcf, 'Name', [prefix ': interaction between ps and lfp'], 'NumberTitle', 'off')
%         
%         if d==1
%             ff = f;
%         end
%         figure(ff);
%         % interaction plot (latency)
%         subplot(2,2,1+(d-1)*2)
%         interaction_plot(para.inter_table_lat(drugtype==d-1, :))
%         xlabel({'pupil size', '(small)'})      
%         ylabel({'pupil size', '(large)'})
%         if d==1
%             title('latency')
%         end
%         % interaction plot (type2reg)
%         subplot(2,2,2+(d-1)*2)
%         interaction_plot(para.type2reg(drugtype==d-1, :), 0)
%         xlabel({'gain change'})      
%         ylabel({'additive change'})
%         if d==1
%             title('stimulus driven LFP')
%         end
%         set(gcf, 'Name', [prefix ': stimulus driven LFP'], 'NumberTitle', 'off')
%         f = f + 3;
%     end
% end