function plot_hmms(hmms, n_state)
% plot results from hmms
%

% input =============
if nargin < 2; n_state = 2; end

% path =======================================
if contains(cd, 'gpfs0')
    mypath = '/gpfs01/nienborg/group';
elseif contains(cd, '172.25.250.112')
    mypath = '\\172.25.250.112\nienborg_group\';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% reference ==================================
load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat'], 'dat');
lenses = length(dat);

% data extraction ========================
lists = zeros(lenses, 2); % is5ht, ismango
varexp = nan(lenses, 2, 3);
nonparamr = nan(lenses, 1);
fr_sort_idx = nan(lenses, n_state, 2);
fr_ses = nan(lenses, n_state, 2);
lent = length(hmms.session(16).cond(1).n_state(n_state).hmm.trnsMat(:));
lene = length(hmms.session(16).cond(1).n_state(n_state).hmm.emtMat(:));
trns = zeros(lenses, lent, 2);
emt = zeros(lenses, lene, 2);
dur = cell(n_state, 2);
sid = cell(n_state, 2);
for i = 1:lenses
    if contains(dat(i).drugname, '5HT')
        lists(i, 1) = 1;
    end
    if contains(dat(i).monkey, 'ma')
        lists(i, 2) = 1;
    end
    for d = 1:2
        % effect size
        nonparamr(i, :) = dat(i).nonparam_ratio;
        
        % variance explained
        for j = 1:3
            varexp(i, d, j) = ...
                hmms.session(i).cond(d).n_state(j).hmm.variance_explained;
        end
        
        % sort by firing rate
        frs = [hmms.session(i).cond(d).n_state(n_state).hmm.state.fr];
        [fr_ses(i, :, d), fr_sort_idx(i, :, d)] = sort(frs);
        
        % transmission, emission matrix
        trns(i, :, d) = hmms.session(i).cond(d).n_state(n_state).hmm.trnsMat(:);
        emt(i, :, d) = hmms.session(i).cond(d).n_state(n_state).hmm.emtMat(:);
        
        % state duration
        for n = 1:n_state
            idx = fr_sort_idx(i, n, d);
            v = length(hmms.session(i).cond(d).n_state(n_state).hmm.state(idx).duration);
            dur{n, d} = [dur{n, d}, v];
            sid{n, d} = [sid{n, d}, i*ones(1, length(v))];
        end            
    end
end

% visualization ========================
close all;
cols = [0 0 0; 1 0 0];
drugnames = {'NaCl', '5HT'};

% variance explained across the number of states
figure;
for k = 1:2
    for d = 1:2
       subplot(2, 2, d + 2*(k-1))
       plot(squeeze(varexp(lists(:,1)==k-1, d, :))', '-', ...
           'color', cols(d,:), 'linewidth', 0.1);
       hold on;
       plot(mean(squeeze(varexp(lists(:,1)==k-1, d, :)), 1), '-', ...
           'color', cols(d,:), 'linewidth', 2) 
       if k==1
           if d==1
               ylabel({'NaCl', 'variance explained'})
               title('base')
           else
               title('drug')
           end
       else
           xlabel('the number of states')
           if d==1
               ylabel({'5HT', 'variance explained'})
           end
       end
       set(gca, 'box', 'off', 'tickdir', 'out')
    end    
end
set(gcf, 'Name', 'variance explained in each state', 'NumberTitle', 'off')

% variance explained & effect size
figure;
for k = 1:2 % NaCl or 5HT
    % plot
    subplot(2, 2, 1 + 2*(k - 1)) 
    unity_scatter(squeeze(varexp(lists(:,1)==k-1, 1, n_state)), ...
        squeeze(varexp(lists(:,1)==k-1, 2, n_state)))
    axis([0 1 0 1])
    if k==1
        ylabel({'NaCl', 'variance explained'})
    else
        xlabel({'baseline', 'variance explained'})
        ylabel({'5HT', 'variance explained'})
    end
    subplot(2, 2, 2 + 2*(k - 1))
    delta = varexp(lists(:,1)==k-1, 1, n_state) - varexp(lists(:,1)==k-1, 2, n_state);
    scatter(delta, nonparamr(lists(:,1)==k-1), 'ok')
    [rr, pp] = corrcoef(delta, nonparamr(lists(:,1)==k-1));
    xx = get(gca, 'XLim');
    yy = get(gca, 'YLim');
    text(xx(1)+0.05*(xx(2)-xx(1)), yy(1)+0.95*(yy(2)-yy(1)), ...
        ['r_{pearson} = ' num2str(rr(1,2))])
    text(xx(1)+0.05*(xx(2)-xx(1)), yy(1)+0.85*(yy(2)-yy(1)), ...
        ['p_{pearson} = ' num2str(pp(1,2))])
    if k==2
        xlabel('\Delta variance explained')
    end
    ylabel('5HT effect size')
    set(gca, 'box', 'off', 'tickdir', 'out')
end
set(gcf, 'Name', 'variance explained', 'NumberTitle', 'off')

% firing rate
try
    figure;
    for k = 1:2 % NaCl or 5HT
        for n = 1:n_state % state
            % plot
            subplot(2, n_state, n + n_state*(k - 1)) 
            unity_scatter(squeeze(fr_ses(lists(:,1)==k-1, n, 1)), ...
                squeeze(fr_ses(lists(:,1)==k-1, n, 2)))
            if k==2
                xlabel('firing rate')
            else
                title(['state = ' num2str(n)])
            end
            if n==1
                ylabel({drugnames{k}, '# sessions'})
            end
        end
    end
    set(gcf, 'Name', 'firing rate', 'NumberTitle', 'off')
catch
   disp('firing rate for the both case is 0') 
end

% transmission matrix
figure;
for k = 1:2 % NaCl or 5HT
    for n = 1:lent % state
        % plot
        subplot(2, lent, n + lent*(k - 1)) 
        unity_scatter(squeeze(trns(lists(:,1)==k-1, n, 1)), ...
            squeeze(trns(lists(:,1)==k-1, n, 2)))
        if k==2
            xlabel('probability)')
        else
            title(['trans param ' num2str(n)])
        end
        if n==1            
            ylabel({drugnames{k}, '# session'})
        end  
        axis([0 1 0 1])
    end
end
set(gcf, 'Name', 'transmission', 'NumberTitle', 'off')

% emission matrix
figure;
for k = 1:2 % NaCl or 5HT
    for n = 1:lene % state
        % plot
        subplot(2, lene, n + lene*(k - 1)) 
        unity_scatter(squeeze(emt(lists(:,1)==k-1, n, 1)), ...
            squeeze(emt(lists(:,1)==k-1, n, 2)))
        if k==2
            xlabel('probability')
        else
            title(['emission param ' num2str(n)])
        end
        if n==1            
            ylabel({drugnames{k}, '# session'})
        end
        axis([0 1 0 1])
    end
end
set(gcf, 'Name', 'emission', 'NumberTitle', 'off')

% duration
figure;
for k = 1:2 % NaCl or 5HT
    for n = 1:n_state % state
        % plot
        subplot(2, n_state, n + n_state*(k - 1)) 
        me = zeros(1, 2);
        v = cell(1, 2);
        for d = 1:2 % base or drug
           v{d} = dur{n, d}(ismember(sid{n, d}, find(lists(:,1)==k-1)));
           h = histogram(v{d}, 'BinWidth', 500);
           h.FaceColor = cols(d, :);
           h.EdgeColor = 'w';
           hold on;
           me(d) = nanmedian(v{d});
        end
        
        % format
        yy = get(gca, 'YLim');
        for d = 1:2
            plot(me(d)*[1 1], yy, '-', 'color', cols(d, :), 'linewidth', 2)
            hold on;
        end
        if k==2
            xlabel('duration (ms)')
        else
            title(['state = ' num2str(n)])
        end
        if n==1            
            ylabel({drugnames{k}, '# transition'})
        end
        set(gca, 'box', 'off', 'tickdir', 'out')   
        
        % stats
        pval = ranksum(v{1}, v{2});
        xx = get(gca, 'XLim');
        text(xx(1)+0.05*(xx(2)-xx(1)), yy(1)+0.95*(yy(2)-yy(1)), ...
            ['p=' num2str(pval)])
    end
end
set(gcf, 'Name', 'duration', 'NumberTitle', 'off')

