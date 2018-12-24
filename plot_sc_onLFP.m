function plot_sc_onLFP(serdata, stmtype)
% test whether a lienar relationship between spike counts and LFP power is changed

close all;
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_useufulfunc'))
addpath(genpath('Z:\Katsuhisa\code\integrated\cbrewer'))

% stimulus type
pairtype = serdata.session(459).results.pairtype;
row = [];
is5ht = [];
for i = 1:length(serdata.session)    
    if serdata.session(i).exist==1 && serdata.session(i).goodunit==1 ...
        && isstruct(serdata.session(i).results) && strcmp(stmtype, serdata.session(i).results.stimulus)  
        row = [row, i];        
        if strcmp(serdata.session(i).results.drugname, '5HT')
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

serdata.session = serdata.session(row);

% pair type
switch pairtype
    case 'drug'
        names = {'baseline', 'drug'};
    case 'sc'
        names = {'low sc', 'high sc'};
    case 'ps'
        names = {'small ps', 'large ps'};
    case 'ps_drug'
        names = {'small ps drug', 'large ps drug'};
end
indnames = {'LFP (uV)', 'delta', 'theta', 'alpha', 'beta', 'gamma'};
lenb = length(indnames);

%%
% extract weights
switch pairtype
    case {'ps', 'sc'}
        is5ht = zeros(1, length(is5ht));
end
ve = cell(1, lenb);
rr = cell(1, lenb);
for b = 1:lenb
    ve{b} = nan(lenr, 2);
    rr{b} = nan(lenr, 2);
    for i = 1:lenr
        for k = 1:2
            X = serdata.session(i).results.cond(k).lfpstm.stm.res(:, 2);
            y = serdata.session(i).results.cond(k).lfpstm.stm.res(:, 2+b);
            % linear regression
            beta = glmfit(X, y, 'normal', 'link', 'identity', 'constant', 'on');
            % variance explained
            pred = glmval(beta, X, 'identity');
            r = corrcoef(y, pred);
            ve{b}(i, k) = r(1, 2)^2;
            % correlation
            r = corrcoef(y, X);
            rr{b}(i, k) = r(1,2);
        end
    end
end

%%
% visualize
fignames = {'variance explained', 'Pearsons r'};
ivals = {ve, rr};
cols = [[0 0 0]; [1 0 0]];
nc = length(unique(is5ht));
for r = 1:2
    figure;
    c = 1;
    pval = nan(1,2);
    for b = 1:lenb
        subplot(2, 3, c)
        med = nan(nc, 2);
            for l = 1:nc
                % scatter
                plot(ivals{r}{b}(is5ht==l-1, 1), ivals{r}{b}(is5ht==l-1, 2), 'o', 'color', cols(l,:), 'markersize', 3)
                hold on;
                
                % stats
                med(l, 1) = median(ivals{r}{b}(is5ht==l-1, 1));
                med(l, 2) = median(ivals{r}{b}(is5ht==l-1, 2));
                pval(l) = signrank(ivals{r}{b}(is5ht==l-1, 1), ivals{r}{b}(is5ht==l-1, 2));
                
                % formatting
                if l==nc
                    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
                    xx = get(gca, 'XLim');
                    yy = get(gca, 'YLim');
                    minima = min([xx, yy]);
                    maxima = max([xx, yy]);
                    plot([minima maxima], [minima maxima], '-', 'color', 0.5*[1 1 1], ...
                        'linewidth', 0.5)
                    axis([minima maxima minima maxima])
                    if minima < 0 && maxima > 0
                        hold on;
                        plot([0 0], [minima maxima], '-', 'color', 0.5*[1 1 1], ...
                            'linewidth', 0.5)
                        hold on;
                        plot([minima maxima], [0 0], '-', 'color', 0.5*[1 1 1], ...
                            'linewidth', 0.5)                        
                        set(gca, 'XTick', [minima 0 maxima], 'YTick', [minima 0 maxima])
                    else
                        set(gca, 'XTick', [minima maxima], 'YTick', [minima maxima])
                    end
                    for j = 1:nc
                        text(minima+0.1*(maxima-minima), maxima-0.1*j*(maxima-minima), ...
                            ['p=' num2str(pval(j)) '(med: ' num2str(med(l, 1)) ' vs ' num2str(med(l, 2)) ')'],...
                            'color', cols(j,:))
                    end
                    % axis
                    title(indnames{b})
                    if c==4
                        xlabel(names{1})
                        ylabel(names{2})
                    end
                end
            end
        c = c + 1;
    end
    set(gcf, 'Name', fignames{r}, 'NumberTitle', 'off')
end