function plot_rcps(rcps)
% plot results from 'pupil_interaction.m' or 'pupil_lfp.m'
% for RC analysis: tuning curve, type II regression
% DRUG, PS, sps: base vs drug, lps: base vs drug

close all;
switch rcps.drugname
    case '5HT'
        dcol = 'r';
    otherwise
        dcol = 'k';
end
fields = rcps.rc.fields;
name = {{'base', 'drug','drug'}, {'lps','sps','ps'}, {3,1,'sps_drug'},{4,2,'lps_drug'}};
for f = 1:length(fields)
    figure(f);
    yall = [];
    for i = 1:4    
        % tuning curve
        mv1 = [];
        mv2 = [];
        if i < 3
            xv1 = [rcps.rc.(['rcsub_' name{i}{1}])(f).results.stm.val];
            yv1 = [rcps.rc.(['rcsub_' name{i}{1}])(f).results.stm.peak];
            xv2 = [rcps.rc.(['rcsub_' name{i}{2}])(f).results.stm.val];
            yv2 = [rcps.rc.(['rcsub_' name{i}{2}])(f).results.stm.peak];
            for k = 1:length(rcps.rc.rcsub_base.results.stm)
                mv1 = [mv1; rcps.rc.(['rcsub_' name{i}{1}])(f).results.stm(k).mean];
                mv2 = [mv2; rcps.rc.(['rcsub_' name{i}{2}])(f).results.stm(k).mean];
            end
        else
            xv1 = [rcps.rc.rcsub(f).each(name{i}{1}).results.stm.val];
            yv1 = [rcps.rc.rcsub(f).each(name{i}{1}).results.stm.peak];
            xv2 = [rcps.rc.rcsub(f).each(name{i}{2}).results.stm.val];
            yv2 = [rcps.rc.rcsub(f).each(name{i}{2}).results.stm.peak];
            for k = 1:length(rcps.rc.rcsub_base.results.stm)
                mv1 = [mv1; rcps.rc.rcsub(f).each(name{i}{1}).results.stm(k).mean];
                mv2 = [mv2; rcps.rc.rcsub(f).each(name{i}{2}).results.stm(k).mean];
            end
        end
        subplot(4,4,i)
        plot(mv1')
        xlabel('time (ms)')
        ylabel('spk/s')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(4,4,i+4)
        plot(mv2')
        xlabel('time (ms)')
        ylabel('spk/s')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(4,4,i+8)
        l = zeros(1,2);
        m = max(yv1);
        yv1 = yv1/m;
        yv2 = yv2/m;
        yall = [yall, yv1, yv2];
        l(1) = plot(xv1, yv1, '-ok');
        hold on;
        l(2) = plot(xv2, yv2, '-or');
        xlabel('orientation (^o)')
        ylabel('total counts per frame')
        if i==2
            legend(l, name{i}{1}, name{i}{2}, 'location', 'northwest')
        else
            legend(l, name{1}{1}, name{1}{2}, 'location', 'northwest')
        end
        legend('boxoff')
        title(name{i}{3})
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % type II regression
        subplot(4,4,i+12)
        scatter(yv1, yv2, 50, 's', ...
            'markerfacecolor', dcol, 'markerfacealpha', 0.4, ...
            'markeredgecolor', 'w', 'markeredgealpha', 0.8)
        hold on;
        minima = min([yv1, yv2])*0.95;
        maxima = max([yv1, yv2])*1.05;
        plot([minima, maxima],[minima, maxima],'--k')
        hold on;
        reg = rcps.rc.type2reg.(name{i}{3});
        plot([minima, maxima], reg(2)*[minima, maxima]+reg(1), ...
            '-', 'color', dcol, 'linewidth',1)
        if i==2
            xlabel(name{i}{1})
            ylabel(name{i}{2})
        else
            xlabel(name{1}{1})
            ylabel(name{1}{2})
        end    
        axis([minima maxima minima maxima])
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
    yrange = [min(yall)*0.95, max(yall)*1.05];
    for i = 1:4
        subplot(4,4,i+8)
        ylim(yrange)
    end
    set(gcf, 'Name', fields{f}, 'NumberTitle', 'off')
end