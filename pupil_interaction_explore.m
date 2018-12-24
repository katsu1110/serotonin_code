function pupil_interaction_explore(anaTps0, anaTps2, varname)
%%
% explore potential interaction of pupil size and drug (5HT)
%
% INPUT: analysis table ... 
%

idx0 = find(contains(anaTps0.varnames, varname));
idx1 = find(contains(anaTps0.varnames, [varname ' base']));
idx2 = find(contains(anaTps0.varnames, [varname ' drug']));

close all;
figure;
animals = {'Kaki', 'Mango', 'Both'};
for a = 1:3 % animals
    subplot(1, 3, a)

    if ~isempty(idx1) && ~isempty(idx2)
        % baseline
        if a < 3
            cond0 = anaTps0.lists(:, 3) == a -1;
        else
            cond0 = ones(size(anaTps0.lists, 1), 1) == 1;
        end
        me = [nanmean(anaTps0.table(cond0, idx1), 1), ...
            nanmean(anaTps0.table(cond0, idx2), 1)];
        sem = [nanstd(anaTps0.table(cond0, idx1), 1), ...
            nanstd(anaTps0.table(cond0, idx2), 1)]./sqrt(sum(cond0==1));
        errorbar([1:2]-0.03, me, sem, 'color', 'k', 'capsize', 0); 
        hold on;

        % 5HT
        if a < 3
            cond2 = anaTps2.lists(:, 2)==1 & anaTps2.lists(:, 3) == a-1;
        else
            cond2 = anaTps2.lists(:, 2)==1;
        end
        me = [nanmean(anaTps2.table(cond2, idx1), 1), ...
            nanmean(anaTps2.table(cond2, idx2), 1)];
        sem = [nanstd(anaTps2.table(cond2, idx1), 1), ...
            nanstd(anaTps2.table(cond2, idx2), 1)]./sqrt(sum(cond2==1));
        errorbar([1:2]+0.03, me, sem, 'color', 'r', 'capsize', 0); 

        % stats
        v = [anaTps0.table(cond0, idx1); anaTps0.table(cond0, idx2); ...
            anaTps2.table(cond2, idx1); anaTps2.table(cond2, idx2)];
        is5ht = [zeros(2*sum(cond0),1); ones(2*sum(cond2), 1)] + 1;
        ps = [zeros(sum(cond0), 1); ones(sum(cond0), 1); ...
            zeros(sum(cond2), 1); ones(sum(cond2), 1)] + 1;
        pvals = anovan(v, {is5ht, ps}, 'model', 'interaction', 'varnames', {'5HT','PS'}, 'display', 'off');

        % format
        xlim([0.8 2.5])
        yy = get(gca, 'YLim');
        yy(2) = yy(2) + 0.5*(yy(2)-yy(1));
        text(0.9, yy(1)+0.95*(yy(2)-yy(1)), ['p(5ht)=' num2str(pvals(1))])
        text(0.9, yy(1)+0.85*(yy(2)-yy(1)), ['p(PS)=' num2str(pvals(2))])
        text(0.9, yy(1)+0.75*(yy(2)-yy(1)), ['p(Intr.)=' num2str(pvals(3))])
        ylim(yy)
        set(gca, 'XTick', 1:2, 'XTickLabel', {'small PS', 'large PS'})
    else
        if a < 3
            cond0 = anaTps0.lists(:, 3) == a - 1;
            cond2 = anaTps2.lists(:, 2)==1 & anaTps2.lists(:, 3) == a-1;
        else
            cond0 = ones(size(anaTps0.lists, 1), 1) == 1;
            cond2 = anaTps2.lists(:, 2)==1;
        end
        me = [nanmean(anaTps0.table(cond0, idx0), 1), ...
            nanmean(anaTps2.table(cond2, idx0), 1)];
        sem = [nanstd(anaTps0.table(cond0, idx0), 1)./sqrt(sum(cond0==1)), ...
            nanstd(anaTps2.table(cond2, idx0), 1)./sqrt(sum(cond2==1))];
        errorbar(1:2, me, sem, 'color', 'k', 'capsize', 0); 
        hold on;

        % stats
        [~, pvals] = ttest2(anaTps0.table(cond0, idx0), anaTps2.table(cond2, idx0));

        % format
        xlim([0.8 2.2])
        yy = get(gca, 'YLim');
        yy(2) = yy(2) + 0.5*(yy(2)-yy(1));
        text(0.9, yy(1)+0.95*(yy(2)-yy(1)), ['p=' num2str(pvals)])
        ylim(yy)
        set(gca, 'XTick', 1:2, 'XTickLabel', {'base', '5HT'})
    end
    set(gca, 'box', 'off', 'tickdir', 'out')
    title(animals{a})
    if a==1
        ylabel(varname)
    end
end
set(gcf, 'Name', [varname '_PS_5HT_Interaction'], 'NumberTitle', 'off')
