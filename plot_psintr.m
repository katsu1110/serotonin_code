function plot_psintr(psintr)
% plot results from 'pupil_interaction.m'

figure;
% pupil time-course
subplot(1,4,1)
t = 1:size(psintr.pupil_timecourse_cntr, 2);
me = mean(psintr.pupil_timecourse_cntr,1);
sem = std(psintr.pupil_timecourse_cntr,[],1)/sqrt(size(psintr.pupil_timecourse_cntr,1));
fill_between(t, me-sem, me+sem,[0 0 0])
hold on;
me = mean(psintr.pupil_timecourse_drug,1);
sem = std(psintr.pupil_timecourse_drug,[],1)/sqrt(size(psintr.pupil_timecourse_drug,1));
fill_between(t, me-sem, me+sem,[1 0 0])
xlabel('time')
ylabel({'z-scored','pupil size'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% interaction table
subplot(1,4,2:3)
imagesc(psintr.inter_table)
text(1-0.25, 1, num2str(psintr.inter_table(1, 1)), 'color', 'r')
text(1-0.25, 2, num2str(psintr.inter_table(2, 1)), 'color', 'r')
text(2-0.25, 1, num2str(psintr.inter_table(1, 2)), 'color', 'r')
text(2-0.25, 2, num2str(psintr.inter_table(2, 2)), 'color', 'r')
title('interaction table (FR)')
set(gca,'XTick',1:2,'XTickLabel',{psintr.drugname,'base'})
set(gca,'YTick',1:2,'YTickLabel',{'s-ps','l-ps'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% model performance
subplot(1,4,4)
bar(1:4, psintr.perf, 'facecolor', 'k', 'edgecolor', 'w')
set(gca,'XTick',1:4,'XTickLabel',{'stm','drug','pupil','interatction'})
xtickangle(45)
ylabel('variance explained')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% rc analysis
if isfield(psintr, 'rc')
    figure;
    subplot(1,3,1)
    rcplot(psintr.rc.rcsub_base.results)
    title('base')
    yy0 = get(gca, 'YLim');
    subplot(1,3,2)
    rcplot(psintr.rc.rcsub_drug.results)
    title(psintr.drugname)
    yy2 = get(gca, 'YLim');
    yy = [min([yy0 yy2]), max([yy0 yy2])];
    subplot(1,3,1)
    ylim(yy)
    subplot(1,3,2)
    ylim(yy)
    
    % type ii regression
    x = [-3 3];
    fieldnames = {'drug', 'ps', 'sps_drug', 'lps_drug'};
    subplot(1,3,3)
    plot([0 0], x, ':k')
    hold on;
    plot(x, [0 0], ':k')
    cols = lines(4);
    p = nan(1, 4);
    for i = 1:4
        hold on;
        p(i) = plot(x, psintr.rc.type2reg.(fieldnames{i})(2)*x ...
            + psintr.rc.type2reg.(fieldnames{i})(1), '-','color',cols(i,:));
    end
    xlim(x)
    title('type II regression')
    legend(p, fieldnames, 'AutoUpdate','off')
    legend('boxoff')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

function rcplot(rc)
wnd = 300;
unique_stm = unique([rc.stm.val]);
nstm = length(rc.stm);
plot([0 wnd+1], 1000*mean(rc.res(:)).*[1 1], ':k')
hold on;
xlim([0 wnd+1])
ylabel({'response', '(mean)'})
xlabel('time (ms)')
set(gca, 'XTick', [0 wnd+1], 'XTickLabel', [0 wnd])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')    
pl = nan(1, nstm);
leg = cell(1, nstm);
col = hsv(nstm);
for n = 1:nstm
    hold on;
    pl(n) = plot(1:wnd, rc.stm(n).mean, '-', 'color', col(n,:));
    leg{n} = ['stm:' num2str(unique_stm(n))];
end
legend(pl, leg, 'AutoUpdate','off')
legend('boxoff')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')