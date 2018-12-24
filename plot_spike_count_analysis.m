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
addpath(genpath([mypath '/Katsuhisa/code/integrated/flexibleModulatedPoisson/code_flexibleModulatedPoisson']))

if nargin < 1; stmtype = {'or', 'co', 'sf', 'sz'}; end
if nargin < 2; analysis = {'all'}; end
if nargin < 3; detrend = 0; end

% load data ====================================
lent = length(stmtype);
% stmdur = 0.45*ones(1, lent);
disp('loading data.....')
for i = 1:lent
    Trmat_temp = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/trialmat/Trmat_' stmtype{i} '.mat']);
    if i==1
        Trmat = Trmat_temp.Trmat;
    else
        Trmat = concatenate_trmat(Trmat, Trmat_temp.Trmat);
    end
%     if strcmp(stmtype{i}, 'rc')
%         stmdur(i) = 2;
%     end
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
                candidate = size(datast{i}.cond(d).lfprel.mat{stms(s)}, 1);
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
                sumat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 1);
                mumat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 2);
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
                candidate = size(datast{i}.cond(d).lfprel.mat{stms(s)}, 1);
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
                susc = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 1);
                musc = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 2);
                
                % detrending
                susc = locdetrend(susc);
                musc = locdetrend(musc);
                
                % noise correlation in each stimulus type
                rr = corrcoef(zscore(susc), zscore(musc));
                nc_temp(s) = 0.5*log((1+rr(1,2))/(1-rr(1,2)));
                
                sumat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 1);
                mumat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 2);
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
                candidate = size(datast{i}.cond(d).lfprel.mat{stms(s)}, 1);
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
                scmat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 1);                
%                 mcmat(:, s) = datast{i}.cond(d).lfprel.mat{stms(s)}(end-mintr+1:end, 2);
            end
            % normal fano factor
            ff(i, 1, d) = nanmean(nanvar(scmat, [], 1)./nanmean(scmat, 1));

            % exponential nonlinearity
            [prshat{i, d}, sff] = fit_flexiblemodulatedpoiss(scmat, 'exp');
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
                me{d} = [me{d}, mean(datast{i}.cond(d).lfprel.mat{stms(s)}(:, 1), 1)];
                % spike count variance
                va{d} = [va{d}, var(datast{i}.cond(d).lfprel.mat{stms(s)}(:, 1), [], 1)];
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