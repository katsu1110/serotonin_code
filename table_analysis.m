function ta = table_analysis(anaT)
%%
% confusioin matrix and dimensionality reduction
%

% inputs
close all;
ndrug = 1;
if strcmp(anaT.pairnames{1,1}, 'base')
    ndrug = 2;
end
lena = length(unique(anaT.lists(:, 3))) + 1;
nvars = size(anaT.table, 2);
animals = {'Kaki', 'Mango', 'Both'};
addpath(genpath('Z:\Katsuhisa\code\integrated\cbrewer\cbrewer'))
cmap = flipud(cbrewer('div', 'RdBu', 101));

% plot
ta.dmvals = nan(size(anaT.lists, 1), 2);
for a = 1:lena
    for k = 1:ndrug
        % condition
        if a < 3
            cond = anaT.lists(:,2)==k-1 & anaT.lists(:, 3)==a-1;
        else
            cond = anaT.lists(:,2)==k-1;
        end
        
        % confusion matrix
        figure(a);
        subplot(1, ndrug, k)
        T = anaT.table(cond, :);
        [ta.animal(a).drug(k).cm, ta.animal(a).drug(k).pm] = ...
            confusion_matrix(T, 'Pearson', 'jet', 0);
        colormap(cmap);
        caxis([-1 1]);
        
        set(gca, 'XTick', 1:nvars, 'XTickLabel', anaT.varnames, 'fontsize', 6)
        set(gca, 'YTick', 1:nvars, 'YTickLabel', anaT.varnames, 'fontsize', 6)
        xtickangle(45); ytickangle(45);
        title(anaT.pairnames{k, 2})
        set(gcf, 'Name', [animals{a}], 'NumberTitle', 'off')
        
        % dimensionality reduction
        figure(a + 100);
        subplot(1, ndrug, k)
        ta.animal(a).drug(k).dm = dimensionality_reduction(T);        
        set(gcf, 'Name', [animals{a}], 'NumberTitle', 'off')
        
        ta.dmvals(cond, :) = ta.animal(a).drug(k).dm.Y(:, 1:2);
    end
end

