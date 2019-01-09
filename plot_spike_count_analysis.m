function plot_spike_count_analysis(stmtype, analysis)
%%
% visualize results from spike count analysis 
%

% path =======================================
if ismember(1, contains(cd, 'gpfs0'))
    mypath = '/gpfs01/nienborg/group';
elseif ismember(1, contains(cd, '\\172.25.250.112'))
    mypath = '//172.25.250.112/nienborg_group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/cbrewer']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/flexibleModulatedPoisson/code_flexibleModulatedPoisson']))

if nargin < 1; stmtype = {'or', 'co', 'sf', 'sz'}; end
if nargin < 2; analysis = {'all'}; end

% load data ====================================
disp('loading data.....')
lent = length(stmtype);
for i = 1:lent
    Trmat_temp = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/trialmat/Trmat_' stmtype{i} '.mat']);
    if i==1
        Trmat = Trmat_temp.Trmat;
    else
        Trmat = concatenate_trmat(Trmat, Trmat_temp.Trmat);
    end
end
disp('data loaded')
pairnames = {'NaCl', '5HT'};

% unit number
lenses = length(Trmat.stmtype);
un = zeros(lenses, 1);
for i = 1:lenses
    low = strfind(Trmat.lfplist{i}{1}, '_');
    if strcmp('ma', Trmat.lfplist{i}{1}(1:2))
        aid = 1000;
    else
        aid = 2000;
    end
    un(i) = aid + str2double(Trmat.lfplist{i}{1}(low(1)+1:low(2)-1));
end

datast = Trmat.LFP_prepro;
lists = [Trmat.goodunit', Trmat.is5ht', Trmat.animal', Trmat.stmtype', un];

close all;
disp(['5HT: ' num2str(sum(lists(:,2)==1)) ' pairs'])
disp(['NaCl: ' num2str(sum(lists(:,2)==0)) ' pairs'])

animals = {'kaki', 'mango', 'both'};
lena = length(animals);

% exclude no good sessions
outs = lists(:, 1)==0;

%     outs(6:9) = 1; % mango
%     outs(1:3) = 1; % mango
%     outs = lists(:,1);
%     outs = 1*(outs < 1);
%     for i = 1:lenses
%         if datast{i}.cond(1).sta.nspk==0 || datast{i}.cond(2).sta.nspk==0
%             outs(i) = 1;
%         end
%     end
datast(outs==1) = [];
lists(outs==1, :) = [];
lenses = size(lists, 1);
disp(['analyzed ses:' num2str(lenses)])

cols = [0 0 0; 1 0 0];
cols2 = [0 0 1; 1 0 0];
j = 0;

%%
% firing rate
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'fr'))==1
    disp('firing rate analysis -------------------------')
    figure(j+1);
    
    % data extraction
    fr = nan(lenses, 2, 2);
    for i = 1:lenses        
        stms = find(datast{i}.stm.vals < 1000);
        lens = length(stms);
        stmdur = 0.45;
        if lists(i, 4)==1
            stmdur = 2;
        end
        for d = 1:2
            % trials x spike counts
            mintr = 10000;
            for s = 1:lens
                candidate = size(datast{i}.cond(d).mat{stms(s)}, 1);
                if candidate < mintr
                    mintr = candidate;
                end
            end
%             if mintr < 5
%                 disp(['session ' num2str(i) ', mintr = ' num2str(mintr)])
%                 continue
%             end
            sumat = nan(mintr, lens);
            mumat = nan(mintr, lens);
            for s = 1:lens
                sumat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 1);
                mumat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 2);
            end
            try
                % firing rate (SU)
                fr(i, 1, d) = nanmean(sumat(:))/stmdur;
                % firing rate (MU)
                fr(i, 2, d) = nanmean(mumat(:))/stmdur;
            catch
                disp(['error in session ' num2str(i)])
                continue
            end
        end
    end
    
    % plot
    clc
    for a = 1:lena
        X = []; Y = []; I = [];
        for d = 1:2
            for k = 1:2
                if a < 3
                    cond = lists(:,3)==a-1 & lists(:,2)==k-1;
                else
                    cond = lists(:,2)==k-1;
                end
                uni = lists(cond, 5);
                x = uniuni(squeeze(fr(cond, d, 1)), uni);
                X = [X; x];
                Y = [Y; uniuni(squeeze(fr(cond, d, 2)), uni)];
                I = [I; zeros(length(x), 1)+(k-1)];
            end
            subplot(lena, 2, d + 2*(a-1))
            unity_scatter(X, Y, I)
            if a==1
                if d==1
                    title('su')
                else
                    title('mu')
                end
            end
        end
        xlabel('FR (base)')
        ylabel('FR (drug)')
    end
    j = j + 1;
end



%%
% GLM analysis ================================================
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'glm'))==1
    disp('GLM analysis -------------------------')
    ndrug = 2;
    paircol = cbrewer('qual', 'Paired', 12);
    varnames = {'stm type', 'ps', 'dps', 'drug', 'drug x ps', 'drug x dps', 'stm pos', ...
        'lfp res', 'csd', 'low-freq', 'gamma', 'label_seq', 'su', 'mu'};
%         mdlnames = {'su', 'mu', 'low-freq'};
%         mdly = [10, 11, 8];
    mdlnames = {'low-freq', 'gamma', 'su'};
    mdly = [10, 11, 12];
    mdlout = {[8:13]};
    lenm = length(mdly);
    lenp = zeros(1, lenm);
    cv = 3;
    cc = cell(1, lenm);
    w = cell(1, lenm);
    lam = 0.084;
    for i = 1:lenses
        % firing rate order
        ss = length(datast{i}.cond(1).mat);
        fr = zeros(1, ss);
        for s = 1:ss
            fr(s) = nanmean(datast{i}.cond(1).mat{s}(:, 1));
        end
        [~, sortidx] = sort(fr);
        X = cell(1, 2);
        ntrall = [0, 0];
        for d = 1:2 % base or drug
            for s = 1:length(datast{i}.cond(1).mat)    % stimulus types
                mat = datast{i}.cond(d).mat{s};
                ntr = size(mat, 1);
                ntrall(d) = ntrall(d) + ntr;
                X{d} = [X{d}; ...
                    [sortidx(s)*ones(ntr, 1), ... % stimulus type
                    mat(:, 3), ... % pupil size
                    mat(:, 4), ... % pupil size derivaticc
                    (d-1)*ones(ntr, 1), ... % base or drug
                    mat(:, 3).*((d-1)*ones(ntr, 1)), ... % interaction
                    mat(:, 4).*((d-1)*ones(ntr, 1)), ... % interaction
                    mat(:, 9), ... % stimulus position in the sequence
                    mat(:, 5), ... % LFP res
                    mat(:, 6), ... % ddLFP
                    mat(:, 7), ... % low-freq
                    mat(:, 8), ... % gamma
                    mat(:, 1), ... % spike counts (SU)
                    mat(:, 2), ... % spike counts (MU)
                    ]]; 
            end
        end

        % k-fold cross-validation
        X = [X{1}; X{2}];
        rng(19891220);
        cvidx1 = crossvalind('Kfold', ntrall(1), cv);
        cvidx2 = crossvalind('Kfold', ntrall(2), cv);
        cvidx = [cvidx1; cvidx2];

        % fit GLM stepwise
        for m = 1:lenm
            if m < 3
                y = 10*log10(X(:, mdly(m)));
            else
                y = log(1+X(:, mdly(m)));
            end
            predictors = zscore(X(:, ~ismember(1:size(X, 2), mdlout{1})));
            lenp(m) = size(predictors, 2);
            r_temp = nan(cv, lenp(m));
            w_temp = nan(cv, lenp(m));
            for v = 1:cv
                for k = 1:lenp(m)
                    % model prediction (stepwise)
                    [B, FitInfo] = lassoglm(predictors(cvidx~=v, 1:k), y(cvidx~=v), 'normal', 'lambda', lam);
                    beta = [FitInfo.Intercept; B];                        
%                     beta = glmfit(predictors(cvidx~=v, 1:k), y(cvidx~=v), 'normal');
                    ypred = glmval(beta, predictors(cvidx==v, 1:k), 'identity');

                    % correlation coefficient
                    rr = corrcoef(y(cvidx==v), ypred);
                    r_temp(v, k) = rr(1,2);
%                         r_temp(v, k) = varexp(y(cvidx==vec(~ismember(vec, v))), ypred);                        
                end
%                     L = [L, FitInfo.LambdaMinDeviance];
                % weight
                w_temp(v, :) = beta(2:end);              
            end 
            [cvscore, cvi] = max(r_temp(:, end));
            cc{m}(i, :) = r_temp(cvi, :);      
            w{m}(i, :) = w_temp(cvi, :);
        end   
    end
    
    % visualize
    c = 1;
    vars = cell(1, m);
    for m = 1:lenm
        for a = 1:lena
            for k = 1:ndrug
                % bar plots ================================================
                figure(j+c);
                if a < 3
                    cond = lists(:,2)==k-1 & lists(:,3)==a-1;
                else
                    cond = lists(:,2)==k-1;
                end
               % individual correlation coefficients
               subplot(lena*ndrug, 2, 1 + 2*(k-1) + ndrug*2*(a-1))
%                    me = nanmean(squeeze(w(cond, :, m)), 1);
%                    sem = nanstd(squeeze(w(cond, :, m)), [], 1)/sqrt(sum(cond));

               plot([0 lenp(m)], [0 0], ':k')
               for p = 1:lenp(m)
                   pval = signrank(w{m}(cond, p));
                   if pval < 0.05/lenp(m)
                        pcol = [1 0 0];
                   else
                       pcol = [0 0 0];
                   end
%                        bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', barcol)
%                        hold on;
%                        errorbar(p, me(p), sem(p), 'color', barcol, 'capsize', 0)
%                        hold on;
                    y = w{m}(cond, p);
                    x = p*ones(size(y));
                    hold on;
                    scatter(x, y, 20, 'o', 'markerfacecolor', pcol, 'markeredgecolor', pcol, ...
                        'markerfacealpha', 0.4, 'markeredgealpha', 0.1)
                    hold on;
                    plot([p-0.25 p+0.25], mean(y)*[1 1], '-', 'color', [0 0.5 0], 'linewidth', 2)
               end
%                    hold on;
               vars{m} = varnames(~ismember(1:lenp(m), mdlout{1}));
%                    violin(squeeze(w(cond, :, m)),'facecolor', pcol, 'edgecolor', [], 'mc', [], 'medc', 'g-')
%                    legend('off')
               if a == 1 && k == 1
                   title(mdlnames{m})
               end
                ylabel({animals{a}, pairnames{k}})
                set(gca, 'XTick', 1:lenp(m), 'XTickLabel', [])
                set(gca, 'box', 'off', 'tickdir', 'out')
               if k == ndrug
                   if a == lena
                       set(gca, 'XTick', 1:lenp(m), 'XTickLabel', vars{m})
                   end
                   xtickangle(45)
               end

               % variance explained
               subplot(lena*ndrug, 2, 2 + 2*(k-1) + ndrug*2*(a-1))
               waterfallchart4mat(cc{m}(cond, :))
               me = nanmean(cc{m}(cond, :), 1);
%                    sem = nanstd(cc{m}(cond, :), [], 1)/sqrt(sum(cond));
               for p = 1:lenp(m)
%                        barcol = paircol(3, :);
                   if p > 1
                       pval = signrank(cc{m}(cond, p-1), cc{m}(cond, p));
                       if pval < 0.05/lenp(m)
                           text(p, me(p)*1.1, '*')
%                                 barcol = paircol(4, :);
                       end
%                            if me(p-1) > me(p) 
%                                bar(p, me(p-1), 'FaceColor', barcol, 'EdgeColor', 'k')
%                                hold on;
%                                bar(p, me(p), 'FaceColor', 'w', 'EdgeColor', 'w')
%                            else
%                                bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', 'k')
%                                hold on;
%                                bar(p, me(p-1), 'FaceColor', 'w', 'EdgeColor', 'w')
%                            end
%                            hold on;
%                            plot([p-1.4 p+0.4], me(p-1)*[1 1], '-k', 'linewidth', 0.25)
                   else
                       pval = signrank(cc{m}(cond, p));
                       if pval < 0.05/lenp(m)
                           text(p, me(p)*1.1, '*')
%                                 barcol = paircol(4, :);
                       end
%                            bar(p, me(p), 'FaceColor', barcol, 'EdgeColor', 'k')
                   end
%                        hold on;
%                        errorbar(p, me(p), sem(p), 'color', barcol, 'capsize', 0)
%                        hold on;
               end
               if a == 1 && k == 1
                   title('correlation coefficient')
               end
               set(gca, 'XTick', 1:lenp(m), 'XTickLabel', [])
               yy = get(gca, 'YLim');
              yy(2) =yy(2)*1.1; 
              ylim([floor(0.7*me(1)) yy(2)])
               set(gca, 'box', 'off', 'tickdir', 'out')
               if k == ndrug
                   if a == lena
                       set(gca, 'XTick', 1:lenp(m), 'XTickLabel', vars{m})
                   end
                   xtickangle(45)
               end                
            end
            set(gcf, 'Name', ['GLM: ' mdlnames{m}], 'NumberTitle', 'off') 

%                 % histogram of weights =============================================
%                figure(j+c+1);
%                if a < 3
%                   cond0 = lists(:,2)==0 & lists(:,3)==a-1;
%                   cond2 = lists(:,2)==1 & lists(:,3)==a-1;
%                else
%                   cond0 = lists(:,2)==0;
%                   cond2 = lists(:,2)==1;
%                end
%                for p = 1:lenp-1
%                    subplot(lena, lenp-1, p + (lenp-1)*(a-1))
%                    [~, pval] = ttest2(squeeze(w(cond0, p, m)), squeeze(w(cond2, p, m)));
%                    histogram(squeeze(w(cond0, p, m)), 'FaceColor', 'k')
%                    hold on;
%                    histogram(squeeze(w(cond2, p, m)), 'FaceColor', 'r')
%                    hold on;
%                    if a == 1
%                        title(varnames{2+p})
%                    end        
%                    if p==1
%                        ylabel(animals{a})
%                    end
%                    set(gca, 'box', 'off', 'tickdir', 'out')
%                    xx = get(gca, 'XLim');
%                    yy = get(gca, 'YLim');
%                    text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.9*(yy(2)-yy(1)), ...
%                        ['p=' num2str(pval)], 'fontsize', 6)
%                    axis([xx yy])
%                end 
        end
%             set(gcf, 'Name', 'GLM weight distribution', 'NumberTitle', 'off') 
        j = j + c + 1;   
    end
end


%%
% noise correlation
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'nc'))==1
    disp('noise corr analysis -------------------------')
    figure(j+1);
    
    % data extraction
    nc = nan(lenses, 2);
    fr = nan(lenses, 2, 2);
    for i = 1:lenses        
        stms = find(datast{i}.stm.vals < 1000);
        lens = length(stms);
        stmdur = 0.45;
        if lists(i, 4)==1
            stmdur = 2;
        end
        for d = 1:2
            % trials x spike counts
            mintr = 10000;
            for s = 1:lens
                candidate = size(datast{i}.cond(d).mat{stms(s)}, 1);
                if candidate < mintr
                    mintr = candidate;
                end
            end
%             if mintr < 5
%                 disp(['session ' num2str(i) ', mintr = ' num2str(mintr)])
%                 continue
%             end
            sumat = nan(mintr, lens);
            mumat = nan(mintr, lens);
            nc_temp = nan(1, lens);
            for s = 1:lens
                % original spike counts
                susc = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 1);
                musc = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 2);
                
                % detrending
                susc = locdetrend(susc);
                musc = locdetrend(musc);
                
                % noise correlation in each stimulus type
                rr = corrcoef(zscore(susc), zscore(musc));
                nc_temp(s) = 0.5*log((1+rr(1,2))/(1-rr(1,2)));
                
                sumat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 1);
                mumat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 2);
            end
            % noise correlation
            nc(i, d) = mean(nc_temp);
            % firing rate (SU)
            fr(i, 1, d) = nanmean(sumat(:))/stmdur;
            % firing rate (MU)
            fr(i, 2, d) = nanmean(mumat(:))/stmdur;
        end
    end
    
    % linear regression for firing rate dependence
    beta = glmfit(abs(squeeze(fr(:, 2, 1) - fr(:, 1, 1))), nc(:, 1));
    for i = 1:lenses
        for d = 1:2
            % firing rate correction
            nc(i, d) = nc(i, d) - glmval(beta, abs(squeeze(fr(i, 2, d) - fr(i, 1, d))), 'identity');
        end
    end
    
    
    % plot
    clc
    for a = 1:lena
        X = []; Y = []; I = [];
        for k = 1:2
            if a < 3
                cond = lists(:,3)==a-1 & lists(:,2)==k-1;
            else
                cond = lists(:,2)==k-1;
            end
            uni = lists(cond, 5);
            x = uniuni(nc(cond, 1), uni);
            X = [X; x];
            Y = [Y; uniuni(nc(cond, 2), uni)];
            I = [I; zeros(length(x), 1)+(k-1)];
        end
        % scatter
        subplot(lena, 2, 1+2*(a-1))
        unity_scatter(X, Y, I)
        xlabel('noise corr (base)')
        ylabel('noise corr (drug)')
        
        % histogram
        subplot(lena, 2, 2+2*(a-1))
        for k = 1:2
            h = histogram(X(I==k-1) - Y(I==k-1));
            h.BinWidth = 0.1;
            h.FaceColor = cols(k, :);
            h.EdgeColor = 'w';
            h.FaceAlpha = 0.4;
            h.EdgeAlpha = 0.4;
            hold on;            
        end
        yy = get(gca, 'YLim');
        for k = 1:2
            % median
            med = median(X(I==k-1) - Y(I==k-1));
            plot(med*[1 1], yy, '-', 'color', cols(k, :))
            hold on;
            
            % stats
            pval = signrank(X(I==k-1), Y(I==k-1));
            text(0.3, yy(1)+(0.95 - 0.1*(k-1))*(yy(2)-yy(1)), ...
                [pairnames{k} ' : p=' num2str(pval)], 'color', cols(k, :), 'fontsize', 6)
        end
        ylim(yy)
        xlabel('\Delta noise corr (base - drug)')
        ylabel('number of units')
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
    j = j + 1;
end

%%
% fano factor
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'ff'))==1
    disp('fano factor analysis -------------------------')
    figure(j+1);
    mdls = {'ff', 'ff_exp'};
    lenm = length(mdls);
    
    % data extraction
    ff = nan(lenses, lenm, 2);
    prshat = cell(lenses, 2);
    for i = 1:lenses        
        stms = find(datast{i}.stm.vals < 1000);
        lens = length(stms);
        for d = 1:2
            % trials x spike counts
            mintr = 10000;
            for s = 1:lens
                candidate = size(datast{i}.cond(d).mat{stms(s)}, 1);
                if candidate < mintr
                    mintr = candidate;
                end
            end
%             if mintr < 5
%                 disp(['session ' num2str(i) ', mintr = ' num2str(mintr)])
%                 continue
%             end
            scmat = nan(mintr, lens);       
%             mcmat = nan(mintr, lens);
            for s = 1:lens
                scmat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 1);                
%                 mcmat(:, s) = datast{i}.cond(d).mat{stms(s)}(end-mintr+1:end, 2);
            end
            % normal fano factor
            ff(i, 1, d) = nanmean(nanvar(scmat, [], 1)./nanmean(scmat, 1));

            % exponential nonlinearity
            try
                [prshat{i, d}, sff] = fit_flexiblemodulatedpoiss(scmat, 'exp');
            catch
                keyboard
            end
            ff(i, 2, d) = mean(sff);
            
%             % noise correlation
%             prshat2 = fit_flexiblemodulatedpoiss(mcmat, 'exp');
%             rr = corrcoef(zscore(prshat{i,d}(1:end-1)), zscore(prshat2(1:end-1)));
%             ff(i, 3, d) = rr(1, 2);            
        end
    end
    
    % plot
    clc
    for a = 1:lena
        for d = 1:lenm % model
            X = []; Y = []; I = [];
            for k = 1:2 % is5HT
                if a < 3
                    cond = lists(:,3)==a-1 & lists(:,2)==k-1;
                else
                    cond = lists(:,2)==k-1;
                end
                uni = lists(cond, 5);
                x = uniuni(squeeze(ff(cond, d, 1)), uni);
                X = [X; x];
                Y = [Y; uniuni(squeeze(ff(cond, d, 2)), uni)];
                I = [I; zeros(length(x), 1)+(k-1)];
            end
            subplot(lena, lenm, d + lenm*(a-1))
            unity_scatter(X, Y, I)
            if a==1
                title(mdls{d})
            end
        end
        xlabel('(base)')
        ylabel('(drug)')
    end
    j = j + 1;
end

%%
% spike count variance vs mean 
if sum(contains(analysis, 'all'))==1 || sum(contains(analysis, 'mevsvar'))==1
    disp('MEvsVAR analysis -------------------------')
    figure(j+1);
    
    % data extraction
    st = cell(1, 2);
    me = cell(1, 2);
    va = cell(1, 2);
    for i = 1:lenses        
        stms = find(datast{i}.stm.vals < 1000);
        lens = length(stms);
        for d = 1:2
            for s = 1:lens
                % stimulus 
                st{d} = [st{d}, i];
                % spike count mean
                me{d} = [me{d}, mean(datast{i}.cond(d).mat{stms(s)}(:, 1), 1)];
                % spike count variance
                va{d} = [va{d}, var(datast{i}.cond(d).mat{stms(s)}(:, 1), [], 1)];
            end
        end
    end
    
    % plot ============================
    for a = 1:lena    
        for k = 1:2
            if a < 3
                cond1 = ismember(st{1}, find(lists(:,2)==k-1 & lists(:,3)==a-1));
                cond2 = ismember(st{2}, find(lists(:,2)==k-1 & lists(:,3)==a-1));
            else
                cond1 = ismember(st{1}, find(lists(:,2)==k-1));
                cond2 = ismember(st{2}, find(lists(:,2)==k-1));
            end
            % baseline
            plot(me{1}(cond1), va{1}(cond1), '.', 'color', cols(1, :))
            hold on;
            % drug
            plot(me{2}(cond2), va{2}(cond2), '.', 'color', cols2(k, :))
            hold on;
        end
        plot([0.01 1000], [0.01 1000], '--', 'color', 0.5*[1 1 1])
        hold on;
    end
    
    % fitting ==========================
    options = optimset('MaxFunEvals',10000,'maxiter',10000);
    x = linspace(0.01, 1000, 100);
    
    % baseline
    c = @(p) cost(p, me{1}, va{1});
    gnoise_base = fminsearch(c, 0, options);
    hold on;
    plot(x, expnonlin(x, gnoise_base), '-k', 'linewidth', 2);
    text(0.05, 500, ['baseline: ' num2str(length(me{1}))], 'fontsize', 6)
    
    % NaCl
    y = me{2}(ismember(st{2}, find(lists(:,2)==0)));
    v = va{2}(ismember(st{2}, find(lists(:,2)==0)));
    c = @(p) cost(p, y, v);
    gnoise_nacl = fminsearch(c, 0, options);
    hold on;
    plot(x, expnonlin(x, gnoise_nacl), '-', 'color', cols2(1,:), 'linewidth', 2);    
    text(0.05, 100, ['NaCl: ' num2str(length(y))], 'color', cols2(1, :), 'fontsize', 6)
    
    % 5HT
    y = me{2}(ismember(st{2}, find(lists(:,2)==1)));
    v = va{2}(ismember(st{2}, find(lists(:,2)==1)));
    c = @(p) cost(p, y, v);
    gnoise_ser = fminsearch(c, 0, options);
    hold on;
    plot(x, expnonlin(x, gnoise_ser), '-', 'color', cols2(2,:), 'linewidth', 2);
    text(0.05, 50, ['5HT: ' num2str(length(y))], 'color', cols2(2, :), 'fontsize', 6)
     
%     if exist('prshat', 'var')
%         x = linspace(0.01, 1000, 100);
%         va = nan(lenses, 100, 2);
%         for i = 1:lenses
%             for d = 1:2
%                 va(i, :, d) = expnonlin(x, prshat{i, d});
%             end
%         end
%         hold on;
%         % baseline
%         plot(x, mean(squeeze(va(:, :, 1)), 1), '-k')
%         hold on;
%         % NaCl
%         plot(x, mean(squeeze(va(lists(:,2)==0, :, 2)), 1), 'color', cols2(1,:))
%         hold on;
%         % 5HT
%         plot(x, mean(squeeze(va(lists(:,2)==1, :, 2)), 1), 'color', cols2(2,:))
%         hold on;
%     end
    
    axis([0.01 1000 0.01 1000])        
    xlabel('mean of spike counts')
    ylabel('variance of spike counts')
    set(gca,'XScale', 'log')
    set(gca,'YScale', 'log')
    set(gca, 'box', 'off', 'tickdir', 'out')
        
    j = j + 1;
end


function Trmat = concatenate_trmat(Trmat1, Trmat2)
fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'LFP_prepro'};
lenf = length(fields);
lenses1 = length(Trmat1.stmtype);
lenses2 = length(Trmat2.stmtype);
for f = 1:lenf
    Trmat1.(fields{f})(lenses1+1:lenses1+lenses2) = Trmat2.(fields{f})(:); 
end
Trmat = Trmat1;

function y = uniuni(x, un)
u = unique(un);
lenu = length(u);
y = nan(lenu, 1);
for i = 1:lenu
    y(i) = mean(x(un==u(i)));
end

function va = expnonlin(x, prshat)
% exponential nonlinearity
va = x + (exp(exp(prshat(end))) - 1)*(x.^2);
% lenv = length(prshat) - 1;
% me = nan(lenv, 1); 
% va = nan(lenv, 1);
% for i = 1:lenv
%     me(i) = exp(prshat(i) + 0.5*prshat(end));
%     va(i) = me(i) + (exp(prshat(end)) - 1)*(me(i)^2);
% end
% [me, i] = sort(me, 'ascend');
% va = va(i);
% me = [range(1); me; range(2)];
% va = [me(1) + (exp(prshat(end)) - 1)*(me(1)^2);
%     va; me(end) + (exp(prshat(end)) - 1)*(me(end)^2)];
% va = interp1(me, va, linspace(me(1), me(end), 100));

function c = cost(p, x, y)
c = mean((y - expnonlin(x, p)).^2);