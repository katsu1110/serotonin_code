function lfps_visualize(lfps, analysis)
%%
% plot results from 'lfps.mat'
%

close all;

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% input ===================
if nargin < 1; end
if nargin < 2; analysis = {'all'}; end

lfpfield = 'LFP_prepro';

% stimulus type ============================
lists = [lfps.goodunit', lfps.is5ht', lfps.animal'];
lfps.cond(1).LFP_prepro(lists(:,1)==0) = [];
lfps.cond(2).LFP_prepro(lists(:,1)==0) = [];
lists(lists(:,1)==0, :) = [];
lenses = size(lists, 1);
switch lfps.stmtype(1)
    case 1 % RC
        stmidx = ones(lenses, 1);
    case 2 % or
        stmidx = zeros(lenses, 2);
        for i = 1:lenses
            [~, maxi] = max(lfps.cond(1).(lfpfield){i}.spk_tu(:,1));
            prfor = lfps.cond(1).(lfpfield){i}.stm.vals(maxi);
            upfor = prfor + 90;
            if upfor > 180
                upfor = upfor - 180;
            end
            [~, mini] = min(abs(lfps.cond(1).(lfpfield){i}.stm.vals - upfor));
            stmidx(i, :) = [mini, maxi];
        end
    case 3 % co
        stmidx = zeros(lenses, 3);
        cos = [0.25, 0.5, 1];
        for i = 1:lenses
            for c = 1:3
                [~, idx] = min(abs(lfps.cond(1).(lfpfield){i}.stm.vals - cos(c)));
                stmidx(i, c) = idx;
            end
        end
    case {4, 5} % sf, sz
        stmidx = zeros(lenses, 2);
        for i = 1:lenses
            [~, maxi] = max(lfps.cond(1).(lfpfield){i}.spk_tu(:,1));
            [~, mini] = min(lfps.cond(1).(lfpfield){i}.spk_tu(:,1));
            stmidx(i, :) = [mini, maxi];
        end
end
ss = size(stmidx, 2);
cols = {[0 0 0; 0 0 1], [0 0 0; 1 0 0]};

% STA ======================================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'sta'))
    % data extraction
    thre = 10;
    vlen = length(lfps.cond(1).(lfpfield){1}.sta.mean);
    mat = {nan(lenses, vlen, ss), nan(lenses, vlen, ss)};
    t = -((vlen-1)/2):((vlen-1)/2);
    st = lfps.cond(1).(lfpfield){end}.sta.t{end};
    st = linspace(t(1), t(end), length(st));
    freq = lfps.cond(1).(lfpfield){end}.sta.f{end};
    spc = {zeros(length(freq), length(st), ss), zeros(length(freq), length(st), ss)};
    nspk = zeros(lenses, 2, ss);
    ntrs = zeros(2, lenses, ss);
    for i = 1:lenses
        for d = 1:2
            nspk(i,d) = lfps.cond(d).(lfpfield){i}.sta.nspk;
            mat{d}(i,:) = lfps.cond(d).(lfpfield){i}.sta.mean;
%             mat{d}(i,:) = lfps.cond(d).(lfpfield){i}.sta.mean ...
%                 - mean(lfps.cond(d).(lfpfield){i}.sta.mean(1:round((vlen-1/2)/2)));
        end
        s1 = lfps.cond(1).(lfpfield){i}.sta.p{end};
        s2 = lfps.cond(2).(lfpfield){i}.sta.p{end};
        if sum(isnan(s1(:)))==0 && sum(isnan(s2(:)))==0
            ntrs(1,i) = lfps.cond(1).(lfpfield){i}.ntr;
            ntrs(2,i) = lfps.cond(2).(lfpfield){i}.ntr;
            spc{lists(i, 2)+1} = spc{lists(i, 2)+1} + (s1 - s2)*sum(ntrs(:,i));
        end
    end

    % visualize
    h = figure;
    clim = zeros(2, 2);
    for k = 1:2
        % condition
        cond = lists(:, 1)==1 & lists(:,2)==k-1 & min(nspk, [], 2) > thre;
        
        % STA ======================================================
        subplot(2, 4, 1 + 4*(k - 1))
%         subplot(1,2,k)
        zeroamp = nan(sum(cond), 2);
        me = zeros(2, length(t)); 
        for d = 1:2
            zeroamp(:, d) = median(mat{d}(cond, t > -5 & t < 5), 2);
            m = mat{d}(cond, :);
            me(d,:) = mean(m, 1);
            sem = std(m, [], 1)/sqrt(sum(cond));
            fill_between(t, me(d,:) - sem, me(d,:) + sem, cols{k}(d, :), 0.4);
            hold on;
            plot(t, me(d,:), '-', 'color', cols{k}(d, :))
            hold on;
        end
        yy = get(gca, 'YLim');
        plot([0 0], yy, ':k')
        ylim(yy)
        xlabel('time from spikes (ms)')
        ylabel('lfp')
        set(gca, 'box', 'off', 'tickdir', 'out')
        
        % stats
        xx = get(gca, 'XLim');
        text(xx(1)+0.05*(xx(2)-xx(1)), yy(1)+0.3*(yy(2)-yy(1)), ['N = ' num2str(sum(cond))])
                
        % zeroamp scatter =====================================
        subplot(2, 4, 2 + 4*(k - 1))
        plot(zeroamp(:,1), zeroamp(:,2), 'ok')
        xx = get(gca, 'XLim');
        yy = get(gca, 'YLim');
        hold on;
        minima = min([xx, yy]);
        maxima = max([xx, yy]);
        plot([minima maxima], [minima maxima], '-', 'color', 0.4*[1 1 1])
        axis([minima maxima minima maxima])
        xlabel('baseline')
        ylabel('drug')
        set(gca, 'box', 'off', 'tickdir', 'out')
        
        % stats
        p = signrank(zeroamp(:,1), zeroamp(:, 2));
        text(minima+0.05*(maxima - minima), minima+0.9*(maxima - minima), 'p(signrank) = ')
        text(minima+0.08*(maxima-minima), minima+0.8*(maxima-minima), num2str(p))
        
        % zeroamp vs firing rate =====================================
        subplot(2, 4, 3 + 4*(k - 1))
        a = zeroamp(:,1)-zeroamp(:,2);
        b = nspk(cond, 1)./sum(ntrs(1, cond)) - nspk(cond, 2)./sum(ntrs(2, cond));
        plot(a, b, 'ok')
        xx = get(gca, 'XLim');
        yy = get(gca, 'YLim');
        xlabel('\Delta zeroamp')
        ylabel('\Delta firing rate')
        set(gca, 'box', 'off', 'tickdir', 'out')
        hold on;
        plot([0 0], yy, ':k')
        hold on;
        plot(xx, [0 0], ':k')
        axis([xx yy])
        hold on;
        beta = glmfit(a, b);
        plot(xx, glmval(beta, xx, 'identity'), '-k')
        
        % stats
        [r, p] = corr(a, b, 'type', 'Spearman');
        text(xx(1)+0.05*(xx(2) - xx(1)), yy(1)+0.8*(yy(2) - yy(1)), ['p(spe.)=' num2str(p)])
        text(xx(1)+0.05*(xx(2) - xx(1)), yy(1)+0.9*(yy(2) - yy(1)), ['r(spe.)=' num2str(r)])
        
        % STA spectrogram ====================================
        subplot(2, 4, 4 + 4*(k - 1))
        imagesc(t, freq, spc{k}/sum(sum(ntrs(:, cond))))
        hold on;
        yy = get(gca, 'YLim');
        plot([0 0], yy, ':w')
        ylim(yy)
        xlabel('time from spikes (ms)')
        ylabel('frequency (Hz)')
        set(gca, 'box', 'off', 'tickdir', 'out')        
        clim(d,:) = caxis;
    end
    clim = [min(clim(:, 1)), max(clim(:, 2))];
    for k = 1:2
       subplot(2, 4, 4 + 4*(k - 1))
       caxis(clim)
       colorbar('eastoutside')
    end
    set(h, 'Name', 'STA', 'NumberTitle', 'off')
end

% Spectrogram ======================================
if ismember(1, contains(analysis, 'all')) || ismember(1, contains(analysis, 'spectrogram'))
    h = figure;
    % data extraction
    window = lfps.cond(1).(lfpfield){end}.window;
    t = lfps.cond(1).(lfpfield){end}.spectrogram.t{end};
    t = linspace(lfps.cond(1).(lfpfield){end}.ts(1), ...
        lfps.cond(1).(lfpfield){end}.ts(end), length(t));
    f = lfps.cond(1).(lfpfield){end}.spectrogram.f{end};
    mat = cell(2,2);
    for i = 1:4
        mat{i} = zeros(length(f), length(t));
    end
    ts = lfps.cond(1).(lfpfield){end}.ts;
    lfpmat = {zeros(lenses, length(ts)), zeros(lenses, length(ts))};
    ntrs = zeros(2, lenses);
    for i = 1:lenses
        for d = 1:2
            if sum(isnan(lfps.cond(d).(lfpfield){i}.spectrogram.p{end}))==0
                m = lfps.cond(d).(lfpfield){i}.spectrogram.p{end}...
                    *lfps.cond(d).(lfpfield){i}.ntr;
                mat{lists(i,2)+1, d} = mat{lists(i,2)+1, d} + m;
                ntrs(d, i) = lfps.cond(d).(lfpfield){i}.ntr;
                lfpmat{d}(i,:) = mean(lfps.cond(d).(lfpfield){i}.lfp{end}(:, 1:length(ts)), 1);
            end
        end
    end
    
    % visualization
    clim = zeros(4,2);
    c = 1;
    for k = 1:2
                
        % LFP trace ====================
        subplot(2, 5, 1+5*(k-1))
        for d = 1:2
            [me, sem] = weighted(lfpmat{d}(lists(:,2)+1==k, :), ntrs(d, lists(:,2)+1==k)');
            fill_between(ts, me - sem, me + sem, cols{k}(d, :), 0.4);
            hold on;
            plot(ts, me, '-', 'color', cols{k}(d,:))
        end
        yy = get(gca, 'YLim');
        hold on;
        plot([0 0], yy, '--k')
        set(gca, 'box', 'off', 'tickdir', 'out')
        xlabel('time after stimulus onset (s)')
        ylabel('lfp')
        
        % spectrogram ==================
        m = cell(1, 2);
        for d = 1:2
            subplot(2, 5, 1+d+5*(k-1))
            m{d} = log10(mat{k, d}/sum(ntrs(d, lists(:,2)+1==k)));
            m{d} = m{d} - repmat(mean(m{d}(:, t < 0), 2), 1, length(t));
            imagesc(t, f, m{d})
            clim(c,:) = caxis;
            c = c + 1;
        end

        % delta ======================
        subplot(2, 5, 4+5*(k-1))
        imagesc(t, f, m{1} - m{2})
        yy = get(gca, 'YLim');
        hold on;
        plot([0 0], yy, '--w')
        title('difference')
        set(gca, 'box', 'off', 'tickdir', 'out')
        
        % power =====================
        subplot(2, 5, 5+5*(k-1))
        for d = 1:2
            y = mean(m{d}(:, t > window{end}(1) & window{end}(2) >= t), 2);
            plot(f, y, '-', 'color', cols{k}(d, :))
            hold on;
        end
        xlabel('frequency (Hz)')
        ylabel('power')
        set(gca, 'box', 'off', 'tickdir', 'out')

    end
    % format
    clim = [min(clim(:,1)), max(clim(:,2))];
    for k = 1:2
        for d = 1:2
            subplot(2, 5, 1+d+5*(k-1))
            caxis(clim)
            yy = get(gca, 'YLim');
            hold on;
            plot([0 0], yy, '--w')
            set(gca, 'box', 'off', 'tickdir', 'out')
            if d==1 && k==1
                ylabel('frequency (Hz)')
            elseif d==1 && k==2
                xlabel('time after stimulus onset (s)')
            end
        end
    end    
    set(h, 'Name', 'Spectrogram', 'NumberTitle', 'off')
end

