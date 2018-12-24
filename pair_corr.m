function pair_corr(anaT, varname)
%%
% ranking of correlation coefficients of one variable with others
%

close all;
figure;
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc'))

% vals
column_num = find(contains(anaT.varnames, varname));
ndrug = 1;
if strcmp(anaT.pairnames{1,1}, 'base')
    ndrug = 2;
end
lena = length(unique(anaT.lists(:, 3))) + 1;
nvars = size(anaT.table, 2);
animals = {'Kaki', 'Mango', 'Both'};
for a = 1:lena
    figure(a);
    for k = 1:ndrug
        % condition
        if a < 3
            cond = anaT.lists(:,2)==k-1 & anaT.lists(:, 3)==a-1;
        else
            cond = anaT.lists(:,2)==k-1;
        end
        
        % plot
        subplot(ndrug, 1, k)
        T = anaT.table(cond, :);
        [cc, pval, idx] = corr_ranking(T, column_num, 'Pearson');
        set(gca, 'XTick', 1:nvars-1, 'XTickLabel', anaT.varnames(idx(2:end)), 'fontsize', 7)
        xtickangle(45);
        disp([animals{a} '; vs ' varname ' -------------------'])
        vars = anaT.varnames(idx);
        for v = 1:nvars - 1
            disp([vars{v+1} '; r=' num2str(cc(v+1)) ', p=' num2str(pval(v+1))])
        end
        if a==1
            title([anaT.pairnames{k, 2} '(n=' num2str(sum(cond==1)) ')'])
        end
    end
    set(gcf, 'Name', [animals{a} '; vs ' varname], 'NumberTitle', 'off')
end