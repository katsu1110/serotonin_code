function plotLFPbyStm(lfps)
% plot results obtained by 'LFPbyStm.m'

close all;

name = {'drug', 'ps_base', 'ps_drug'};
winlab = {'baseline (-0.2 - 0s)', 'stimulus evoked (0.05 - 0.3s)', 'sustained (0.35 - 2s)'};
leg = {{'baseline', lfps.drugname}, {'small ps (base)', 'large ps (base)'}, ...
    {['small ps (' lfps.drugname ')'], ['large ps (' lfps.drugname ')']}};
wnd = lfps.drug.cond(1).lfpstm.wnd;
wint = -wnd:0.001:wnd;
for a = 1:3
    % LFP traces
    h = figure;
    l = zeros(1,2);
    me = lfps.(name{a}).cond(1).lfpstm.lfp_stm.mean(end, :);
    sem = lfps.(name{a}).cond(1).lfpstm.lfp_stm.sem(end, :);
    fill_between(lfps.(name{a}).cond(1).lfpstm.ts, me-sem, me+sem,...
        [0 0 0], 0.4)
    hold on;
    l(1) = plot(lfps.(name{a}).cond(1).lfpstm.ts, me, '-k');
    hold on;
    me = lfps.(name{a}).cond(2).lfpstm.lfp_stm.mean(end, :);
    sem = lfps.(name{a}).cond(2).lfpstm.lfp_stm.sem(end, :);
    fill_between(lfps.(name{a}).cond(2).lfpstm.ts, me-sem, me+sem,...
        [1 0 0], 0.4)
    hold on;
    l(2) = plot(lfps.(name{a}).cond(2).lfpstm.ts,me, '-r');
    hold on;
    yy = get(gca, 'YLim');
    plot([0 0], yy, '--k')
    xlim([lfps.(name{a}).cond(2).lfpstm.ts(1), lfps.(name{a}).cond(2).lfpstm.ts(end)])
    ylim(yy)
    legend(l, [leg{a}{1} ':' num2str(lfps.(name{a}).cond(1).lfpstm.ntr)], ...
        [leg{a}{2} ':' num2str(lfps.(name{a}).cond(2).lfpstm.ntr)], 'location', 'southeast')
    legend('boxoff')
    xlabel('time (sec)')
    ylabel('LFP (uV)')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': LFP traces'], 'NumberTitle', 'off')
    
    % power of LFP traces
    h = figure;
    l = zeros(1,2);
    for u = 1:3
        subplot(1,3,u)
        l(1) = plot(lfps.(name{a}).cond(1).lfpstm.period(u).freq(end,:), ...
            lfps.(name{a}).cond(1).lfpstm.period(u).pow_avg(end,:), '-k');
        hold on;
        l(2) = plot(lfps.(name{a}).cond(2).lfpstm.period(u).freq(end,:), ...
            lfps.(name{a}).cond(2).lfpstm.period(u).pow_avg(end,:), '-r');
        hold on;
        xlim([lfps.(name{a}).cond(2).lfpstm.period(u).freq(end,1)...
            lfps.(name{a}).cond(2).lfpstm.period(u).freq(end,end)])
        xlabel('frequency (Hz)')
        ylabel('power (log)')
        title(winlab{u})
        set(gca,'YScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
    end   
    legend(l, [leg{a}{1} ':' num2str(lfps.(name{a}).cond(1).lfpstm.ntr)], ...
        [leg{a}{2} ':' num2str(lfps.(name{a}).cond(2).lfpstm.ntr)], 'location', 'southeast')
    legend('boxoff')
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': LFP power'], 'NumberTitle', 'off')
    
    % spike-triggered LFP
    h = figure;
    l = zeros(1,2);
    for u = 1:3
        subplot(1,3,u)
        me = lfps.(name{a}).cond(1).lfpstm.period(u).stlfp.avg_stlfp;
        sem = lfps.(name{a}).cond(1).lfpstm.period(u).stlfp.sem_stlfp;
        fill_between(wint, me-sem, me+sem, zeros(1,3), 0.4)
        hold on;
        l(1) = plot(wint, me, '-k');
        hold on;
        me = lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.avg_stlfp;
        sem = lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.sem_stlfp;
        fill_between(wint, me-sem, me+sem, [1,0,0], 0.4)
        hold on;
        l(2) = plot(wint, me, '-r');
        hold on;
        yy = get(gca, 'YLim');
        plot([0 0], yy, '--k')
        xlim([-wnd wnd])
        ylim(yy)
        xlabel('time (sec)')
        ylabel('LFP (uV)')
        title(winlab{u})
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
    end   
    legend(l, [leg{a}{1} ':' num2str(lfps.(name{a}).cond(1).lfpstm.ntr)], ...
        [leg{a}{2} ':' num2str(lfps.(name{a}).cond(2).lfpstm.ntr)], 'location', 'southeast')
    legend('boxoff')
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': spike-triggered average LFP'], 'NumberTitle', 'off')
        
    % power of stLFP
    h = figure;
    l = zeros(1,2);
    for u = 1:3
        subplot(1,3,u)
        l(1) = plot(lfps.(name{a}).cond(1).lfpstm.period(u).stlfp.FREQ(:,end), ...
            lfps.(name{a}).cond(1).lfpstm.period(u).stlfp.POW(:,end), '-k');
        hold on;
        l(2) = plot(lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.FREQ(:,end), ...
            lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.POW(:,end), '-r');
        hold on;
%         xlim([lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.FREQ(1,end)...
%             lfps.(name{a}).cond(2).lfpstm.period(u).stlfp.POW(end,end)])
        xlabel('frequency (Hz)')
        ylabel('power (log)')
        title(winlab{u})
        set(gca,'YScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
    end   
    legend(l, [leg{a}{1} ':' num2str(lfps.(name{a}).cond(1).lfpstm.ntr)], ...
        [leg{a}{2} ':' num2str(lfps.(name{a}).cond(2).lfpstm.ntr)], 'location', 'southeast')
    legend('boxoff')
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': stLFP power'], 'NumberTitle', 'off')
     
    % Spectrogram 
    h = figure;
    crange = [0, 1];
    for u = 1:3
        for k = 1:2                
                subplot(3, 3, u+3*(k-1))
                imagesc(lfps.(name{a}).cond(k).lfpstm.period(u).spectrogram.t{end}, ...
                    lfps.(name{a}).cond(k).lfpstm.period(u).spectrogram.f{end},...
                    mag2db(lfps.(name{a}).cond(k).lfpstm.period(u).spectrogram.S{end})')
                axis xy;
                colormap('jet')
                c = caxis;
                colorbar('off')            
                if crange(1) > c(1)
                    crange(1) = c(1);
                end
                if crange(2) < c(2)
                    crange(2) = c(2);
                end           
                caxis(crange)
                if u==1 && k==2
                    message = sprintf([leg{a}{2} ' \n frequency (Hz)']);
                    ylabel(message);
                else
                    ylabel('')
                end
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end

        % difference
        crange = [-1 1];
        subplot(3, 3, u + 6)
        imagesc(lfps.(name{a}).cond(1).lfpstm.period(u).spectrogram.t{end}, ...
                lfps.(name{a}).cond(2).lfpstm.period(u).spectrogram.f{end},...
                (mag2db(lfps.(name{a}).cond(1).lfpstm.period(u).spectrogram.S{end})...
                - mag2db(lfps.(name{a}).cond(2).lfpstm.period(u).spectrogram.S{end}))')
        axis xy;
        colormap('jet')
        colorbar('off')

        c = caxis;
        if crange(1) > c(1)
            crange(1) = c(1);
        end
        if crange(2) < c(2)
            crange(2) = c(2);
        end
        caxis(crange)
        title('')
        xlabel('time [s]');
        message = sprintf([leg{a}{1} ' - ' leg{a}{2} ' \n frequency (Hz)']);
        ylabel(message);
        xlabel('')
        ylabel('')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': spectrogram'], 'NumberTitle', 'off')
    
     % Spike-LFP coherency ===================== 
    h = figure;
    for u = 1:3
        subplot(2, 3, u)
        plot(lfps.(name{a}).cond(1).lfpstm.period(u).coherence.f{end}, ...
            lfps.(name{a}).cond(1).lfpstm.period(u).coherence.C{end}, '-k')
        hold on;
        plot(lfps.(name{a}).cond(2).lfpstm.period(u).coherence.f{end}, ...
            lfps.(name{a}).cond(2).lfpstm.period(u).coherence.C{end}, '-r')
        ylabel('coherence');
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        subplot(2, 3, u+3)
        plot(lfps.(name{a}).cond(1).lfpstm.period(u).coherence.f{end}, ...
            zeros(1, length(lfps.(name{a}).cond(1).lfpstm.period(u).coherence.f{end})), '--k')
        hold on;
        plot(lfps.(name{a}).cond(1).lfpstm.period(u).coherence.f{end}, ...
            lfps.(name{a}).cond(1).lfpstm.period(u).coherence.C{end} - ...
            lfps.(name{a}).cond(2).lfpstm.period(u).coherence.C{end}, '-b')
        xlabel('frequency (Hz)');
        ylabel('\Delta coherence')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
    set(h, 'Name', [lfps.monkey num2str(lfps.id) '; ' name{a} ': coherence'], 'NumberTitle', 'off')
    
end