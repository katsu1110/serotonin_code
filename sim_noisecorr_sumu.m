function [nc_su, nc_mu] = sim_noisecorr_sumu(varargin)
%%
% simulation for noise correlation between single unit and multi unit
%

% default input arguments
n_pop = 10;
mfr = randi(50, 1, n_pop);
n_tr = 100;
sig = 0;
plot_flag = 0;
j = 1;              
while  j <= length(varargin)
    switch varargin{j}
        case 'ntr'
            n_tr = varargin{j+1};
            j = j + 2;
        case 'nneuron'
            n_pop = varargin{j+1};
            j = j + 2;
        case 'sig'
            sig = varargin{j+1};               
             j = j + 2;
        case 'plot'
            plot_flag = 1;
            j = j + 1;
    end
end

% trial-by-trial spikecounts
trmat = nan(n_tr, n_pop);
for i = 1:n_pop
    trmat(:, i) = poissrnd(mfr(i), n_tr, 1);
end

% shared noise
sn = normrnd(0, sig, n_tr, 1);
trmat = round(trmat + repmat(sn, 1, n_pop));

% noise correlation 
nc_su = nan(n_pop, n_pop); % between single units
nc_mu = nan(n_pop, 1); % between multi units
for i = 1:n_pop
    % vs other single units
    for k = 1:n_pop
        rr = corrcoef(zscore(trmat(:, i)), zscore(trmat(:, k)));
        nc_su(i, k) = 0.5*log((1 + rr(1, 2))/(1 - rr(1, 2)));
    end
    % vs multi units
    mu = sum(trmat(:, ~ismember(1:n_pop, i)), 2);
    rr = corrcoef(zscore(trmat(:, i)), zscore(mu));
    nc_mu(i, 1) = 0.5*log((1 + rr(1, 2))/(1 - rr(1, 2)));
end

% visualize
if plot_flag
    figure;
    subplot(1,2,1)
    imagesc(nc_su)
    xlabel('neurons')
    ylabel('neurons')
    colormap(hot)
    c = colorbar('southoutside');
    c.Label.String = 'noise correlation (vs su)';

    subplot(1,2,2)
    for i = 1:n_pop
        plot([1 2], [nc_mu(i), mean(nc_su(i, ~ismember(1:n_pop, i)), 2)], '-ok')
        hold on;
    end
    xlim([0.75 2.25])
    set(gca, 'XTick', [1 2], 'XTickLabel', {'vs mu', 'mean(vs su)'})
    ylabel('noise correlation (vs mu)')
end
