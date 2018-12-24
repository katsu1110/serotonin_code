function irasa_visualize(stmtype, splittype, internal)
%%
% visualize results from IRASA decomposition analysis
%

% input ====================
if nargin < 1; stmtype = 'rc'; end
if nargin < 2; splittype = 'drug'; end
if nargin < 3; internal = 0; end

% path ==========================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/cbrewer/cbrewer']))

% load dataset ====================
ndrug = 1;
switch splittype
    case 'drug'
        ider = '';
        ndrug = 2;
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        ider = '_sc';
        pairnames = {'small sc', 'large sc'};
    case 'pupil'
        ider = '_ps';
        pairnames = {'small ps', 'large ps'};
end
if internal==0
    Irasas = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas' ider '/Irasas_' stmtype '.mat'],...
        'Irasas');
    Irasas = Irasas.Irasas;
else
    Irasas = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas' ider '/iIrasas_' stmtype '.mat'],...
        'iIrasas');
    Irasas = Irasas.iIrasas;
end

disp(['stmtype: ' stmtype])
disp(['mango: ' num2str(sum(Irasas.animal==1))...
    ', kaki: ' num2str(sum(Irasas.animal==0))])
disp(['total ' num2str(length(Irasas.goodunit)) ...
    ', good unit: ' num2str(sum(Irasas.goodunit==1))])

% is5ht = Irasas.is5ht(Irasas.goodunit==1);
% ismango = Irasas.animal(Irasas.goodunit==1);
% Irasas.pairs(Irasas.goodunit==0) = [];
is5ht = Irasas.is5ht;
ismango = Irasas.animal;
lenses = length(Irasas.pairs);
xrange = [3 48];
if ndrug==1
    is5ht = zeros(1, length(is5ht));
end

% stimulus type ============================
lasti = 3;
switch stmtype
    case 'rc'
        stmidx = ones(lenses, 1);
        stmlab = {'ORxRC'};
    case 'or'
        stmlab = {'anti-pref', 'pref'};
        stmidx = zeros(lenses, 2);
        for i = 1:lenses      
            try
                [~, maxi] = max(Irasas.pairs{i}.enc{1}{1}.mean);
                [~, mini] = min(Irasas.pairs{i}.enc{1}{1}.mean);
            catch
                [~, maxi] = max(Irasas.pairs{i}.enc{1}{end}.mean);
                [~, mini] = min(Irasas.pairs{i}.enc{1}{end}.mean);
            end
            stmidx(i, :) = [mini, maxi];
        end
    case 'co'
        stmidx = zeros(lenses, 3);
        stmlab = {'0.25', '0.5', '1'};
        cos = [0.25, 0.5, 1];
        for i = 1:lenses
            unistm = Irasas.pairs{i}.enc{1}{1}.unistm;
            for c = 1:3
                [~, idx] = min(abs(unistm - cos(c)));
                stmidx(i, c) = idx;
            end
        end
    case {'sf', 'sz'}
            stmidx = zeros(lenses, 2);
            stmlab = {'min res stm', 'max res stm'};
            for i = 1:lenses
                [~, maxi] = max(Irasas.pairs{i}.enc{1}{1}.mean);
                [~, mini] = min(Irasas.pairs{i}.enc{1}{1}.mean);
                stmidx(i, :) = [mini, maxi];
            end    
end
ss = size(stmidx, 2);

% data extraction =======================
bandnames = {'low-freq', 'alpha', 'beta', 'gamma'};
bandrange = {[3, 10], [8, 12], [15, 25], [40 48]};
lenb = length(bandnames);
paranames = {'fr', 'beta3-10', 'cons3-10', 'beta15-25', 'cons15-48', ...
    'frac 3-7', 'frac 8-12', 'frac 15-25', 'frac 40-48', ...
    'osci 3-7', 'osci 8-12', 'osci 15-25', 'osci 40-48'};
statsrange = {[3, 10], [11, 24], [25, 48]};
lensr = length(statsrange);
lenp = length(paranames);
encnames = {'mean', 'reliability', 'selectivity', 'snr2', 'discriminability', 'entropy', 'c-entropy', 'MI'};
lene = length(encnames);
ts = Irasas.pairs{end}.ts;
freq = Irasas.pairs{end}.freq';
lent = length(ts);
lenf = length(freq);
fields = {'frac', 'osci'};
for f = 1:2 % fractal or oscillation
    for d = 1:2 % base or drug
        for s = 1:ss % stimulus types
            para.cond(d).stm(s).frac = zeros(lenses, lenf, lent);
            para.cond(d).stm(s).osci = zeros(lenses, lenf, lent);
            para.cond(d).stm(s).frac_v = zeros(lenses, lenf);
            para.cond(d).stm(s).osci_v = zeros(lenses, lenf);
            para.cond(d).stm(s).mix_v = zeros(lenses, lenf);
        end
        para.cond(d).enc = nan(lenses, lenp, lene);
    end
end
for i = 1:lenses % sessions
    for d = 1:2 % base or drug
        for s = 1:ss % stimulus types
            % fractal
            para.cond(d).stm(s).frac(i, :, :) = ...
                squeeze(Irasas.pairs{i}.frac{d}(:, :, stmidx(i, s)));
            para.cond(d).stm(s).frac_v(i, :) = squeeze(nanmean(para.cond(d).stm(s).frac(i, :, end-lasti+1:end), 3));
%             para.cond(d).stm(s).frac_r(i, :) = para.cond(d).stm(s).frac_v(i, :) - squeeze(nanmean(...
%                 para.cond(d).stm(s).frac(i, :, ts <= 0), 3));
            
            % oscillation
            para.cond(d).stm(s).osci(i, :, :) = ...
                squeeze(Irasas.pairs{i}.osci{d}(:, :, stmidx(i, s)));
            para.cond(d).stm(s).osci_v(i, :) = squeeze(nanmean(para.cond(d).stm(s).osci(i, :, end-lasti+1:end), 3));
%             para.cond(d).stm(s).osci_r(i, :) = para.cond(d).stm(s).osci_v(i, :) - squeeze(nanmean(...
%                 para.cond(d).stm(s).osci(i, :, ts <= 0), 3));

            % mix
            para.cond(d).stm(s).mix_v(i, :) = squeeze(nanmean(Irasas.pairs{i}.frac{d}(:, end-lasti+1:end, stmidx(i, s)), 2));
        end
        
        % encoding analysis
        if ~strcmp(stmtype, 'rc')
            for j = 1:lenp
                for e = 1:lene
                    try
                        if strcmp('mean', encnames{e}) || strcmp('snr2', encnames{e})
                            temp = mean(Irasas.pairs{i}.enc{d}{j}.(encnames{e}));
                        elseif strcmp('selectivity', encnames{e}) && strcmp(stmtype, 'or')
                             temp = Irasas.pairs{i}.enc{d}{j}.unique.circularvariance;
                        elseif  strcmp('entropy', encnames{e})
                            temp = Irasas.pairs{i}.enc{d}{j}.metabcost(1);
                        elseif  strcmp('c-entropy', encnames{e})
                            temp = Irasas.pairs{i}.enc{d}{j}.metabcost(2);
                        elseif  strcmp('MI', encnames{e})
                            temp = Irasas.pairs{i}.enc{d}{j}.metabcost(3);
                        else
                            temp = Irasas.pairs{i}.enc{d}{j}.(encnames{e});
                        end
                        para.cond(d).enc(i, j, e) = temp;
                    catch
                        para.cond(d).enc(i, j, e) = nan;
                    end
                end
            end
        else
            for j = 1:lenp
                para.cond(d).enc(i, j, 1) = nanmedian(Irasas.pairs{i}.tumats{d}(:, 1+j), 1);
            end
        end
    end
end

% visualize =======================
% average fractal & oscillation
animal_lab = {'kaki', 'mango', 'both'};
close all;
fign = 1;
% for r = 1:2
%     for f = 1:2
%         figure(fign);
%         clim = {[0 0], [0 0]};
%         for s = 1:ss % stimulus type
%             % time vs freq
%             for a = 1:2 % animal
%                 m = cell(2, 2); 
%                 for c = 1:2 % is5ht
%                     for d = 1:2 % base & drug
%                         subplot(2*ss*3, 3, d+3*(s-1)+ss*3*(c-1)+2*ss*(a-1))
%                         if a < 3
%                             cond = is5ht==c-1 & ismango==a-1;
%                         else
%                             cond = is5ht==c-1;
%                         end
%                         m{c, d} = squeeze(nanmean(para.cond(d).stm(s).(fields{f})(cond, :, :), 1));
%                         if f==1
%                             m{c, d} = log(m{c, d});
%                         end
%                         if r==2
%                             m{c, d} = m{c, d} - repmat(mean(m{c, d}(:, ts <= 0), 2), 1, length(ts));
%                         end
%                         imagesc(ts, freq, m{c, d})
%                         colormap(jet)
%                         c_temp = caxis;
%                         if c_temp(1) < clim{1}(1)
%                             clim{1}(1) = c_temp(1);
%                         end
%                         if c_temp(2) > clim{1}(2)
%                             clim{1}(2) = c_temp(2);
%                         end
%                         cb = colorbar('northoutside');
%                         cb.Label.String = tlab{d};
%                         hold on;
%                         yy = get(gca, 'YLim');
%                         plot([0 0], yy, '--w')
%                         set(gca, 'box', 'off', 'tickdir', 'out')
%                         if d==1
%                             ylabel({animal_lab{a}, pairs{c}, stmlab{s}})
%                             if c==2
%                                 xlabel('time after stimulus onset (sec)')
%                             end
%                         end
%                     end
%                     
%                     % difference
%                     subplot(2*ss*3, 3, 3+3*(s-1)+ss*3*(c-1)+2*ss*(a-1))
%                     imagesc(ts, freq, m{c, 1} - m{c, 2})
%                     c_temp = caxis;
%                     if c_temp(1) < clim{2}(1)
%                         clim{2}(1) = c_temp(1);
%                     end
%                     if c_temp(2) > clim{2}(2)
%                         clim{2}(2) = c_temp(2);
%                     end
%                     cb = colorbar('northoutside');
%                     cb.Label.String = '\Delta';
%                     hold on;
%                     yy = get(gca, 'YLim');
%                     plot([0 0], yy, '--w')
%                     set(gca, 'box', 'off', 'tickdir', 'out')      
%                     
%                     % fix color range
%                     for d = 1:2
%                         subplot(2*ss*3, 3, d+3*(s-1)+ss*3*(c-1)+2*ss*(a-1))
%                         caxis(clim{1})                
%                     end
%                     subplot(2*ss*3, 3, 3+3*(s-1)+ss*3*(c-1)+2*ss*(a-1))
%                     caxis(clim{2}) 
%                 end
%             end
%             fign = fign + 1;
%             if r==1
%                 set(gcf, 'Name', [stmtype '_' fields{f}], 'NumberTitle', 'off')
%             else
%                 set(gcf, 'Name', [stmtype '_' fields{f} '_response'], 'NumberTitle', 'off')
%             end
%         end
%     end
% end

% batch
nrow = floor(sqrt(lenses));
ncol = ceil(lenses/nrow);
bcols = [0 0 0; 0 1 1; 1 0 0; 1 0 1];
for s = 1:ss
    for i = 1:lenses
        for d = 1:2
            % mixed signal
            figure(fign);
            subplot(nrow, ncol, i)
            hold on;
            plot(freq, para.cond(d).stm(s).mix_v(i, :), '-', 'color', bcols(1+2*(d-1), :))
            hold on;
            plot(freq, para.cond(d).stm(s).frac_v(i, :), '-', 'color', bcols(1+2*(d-1), :))
            
            % oscillatory signal
            figure(fign + 1);
            subplot(nrow, ncol, i)
            hold on;
            plot(freq, para.cond(d).stm(s).osci_v(i, :), '-', 'color', bcols(1+2*(d-1), :))
        end
        for d = 1:2
            figure(fign -1 + d)
            xlim(xrange)
            set(gca, 'box', 'off', 'tickdir', 'out')
            if d==1
                set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
            else
                set(gca, 'XScale', 'log');
            end
        end
    end
    figure(fign);
    set(gcf, 'Name', 'irasa (fractal + oscillation)', 'NumberTitle', 'off')
    savefig(gcf, [mypath ...
        '/Katsuhisa/serotonin_project/LFP_project/Manuscript/Figures/Fig3_IRASA/raw_figs/' ...
        'mix_individuals_' splittype '.fig'])
    figure(fign+1);
    set(gcf, 'Name', 'irasa (oscillation)', 'NumberTitle', 'off')
    savefig(gcf, [mypath ...
        '/Katsuhisa/serotonin_project/LFP_project/Manuscript/Figures/Fig3_IRASA/raw_figs/' ...
        'osci_individuals_' splittype '.fig'])
    fign = fign + 2;   
end

% last bins
% cols = cbrewer('qual', 'Set1', 9);
% if ss == 3
%     lcols = [cols(2,:); cols(3,:); cols(1,:)];
% elseif ss == 1
%     lcols = [0 0 0; 1 0 0];
% else
%     lcols = [cols(2,:); cols(1,:)];
% end

% lspc = {'-', '--'};
lcols = [0 0 0; 1 0 0];
lspc = {'-', '-'};
% post = {'v', 'r'};
for f = 1:2
    figure(fign);
    for s = 1:ss % stimulus type
        for a = 1:3 % animal
            pairs = cell(lenb, 2);
            for c = 1:ndrug % is5ht                
                if a < 3
                    cond = is5ht==c-1 & ismango==a-1;
                else
                    cond = is5ht==c-1;
                end
                
%                 % stats
%                 pval = ones(1, lensr);
%                 for r = 1:lensr
%                     frange = statsrange{r}(1) <= freq & freq <= statsrange{r}(2);
%                     [~, pval(r)] = ttest(nanmean(para.cond(1).stm(s).([fields{f} '_v'])(cond, frange), 1), ...
%                         nanmean(para.cond(2).stm(s).([fields{f} '_v'])(cond, frange), 1));
%                 end
                
                % overall power
                subplot(3*ndrug, 2, 1+2*(c-1)+ndrug*2*(a-1))
                me = cell(1, 2); 
                for d = 1:2 % base or drug
                    me{d} = nanmean(para.cond(d).stm(s).([fields{f} '_v'])(cond, :), 1);
                    sem = nanstd(para.cond(d).stm(s).([fields{f} '_v'])(cond, :), [], 1)...
                        /sqrt(sum(is5ht==c-1));
                    fill_between(freq, me{d} - sem, me{d} + sem, lcols(d, :), 0.1)
                    hold on;
                    plot(freq, me{d}, lspc{d}, 'color', lcols(d,:))
                    hold on;
%                     if f==1
%                         set(gca, 'XScale', 'log')
%                     end
                end
                set(gca, 'box', 'off', 'tickdir', 'out')
                xlim(xrange)
                ylabel({animal_lab{a}, fields{f}, pairnames{c, 2}})
                if c==1 && a==1
                    title('analysis window')
                elseif a==3 && c==2
                    xlabel('frequency (Hz)')
                end                
                
                % difference 
                subplot(3*ndrug, 2, 2+2*(c-1)+ndrug*2*(a-1))
                plot(freq, me{1}-me{2}, '-', 'color', [0 1 0]/2)
%                 if f==1
%                     set(gca, 'XScale', 'log')
%                 end
                xlim(xrange)
                set(gca, 'box', 'off', 'tickdir', 'out')
                if a==1 && c==1
                    title('difference')
                end
                
%                 % stats         
%                 for b = 1:length(freq)
%                    xf = nanmean(para.cond(1).stm(s).([fields{f} '_v'])(cond, b), 2);
%                    yf = nanmean(para.cond(2).stm(s).([fields{f} '_v'])(cond, b), 2);
%                    nans = isnan(xf) | isnan(yf);
%                    xf(nans) = []; yf(nans) = [];
%                    pairs{b, c} = [xf, yf];
%                 end
                for b = 1:lenb
                   frange = bandrange{b}(1) <= freq & freq <= bandrange{b}(2); 
                   xf = nanmean(para.cond(1).stm(s).([fields{f} '_v'])(cond, frange), 2);
                   yf = nanmean(para.cond(2).stm(s).([fields{f} '_v'])(cond, frange), 2);
                   nans = isnan(xf) | isnan(yf);
                   xf(nans) = []; yf(nans) = [];
                   pairs{b, c} = [xf, yf];
                end         
            end

            % stats         
            sfreq = freq(freq >= 3 & freq <= 48);
%             start = freq(1)*[1 1];
            for b = 1:lenb
                % pair tests
                stats = pair_tests(pairs{b, 1}, pairs{b, 2});

                % condition 1
                subplot(3*ndrug, 2, 1+ndrug*2*(a-1))
                yl = get(gca, 'YLim');
                disp(['cond 1:' bandnames{b} ', ' num2str(a) ', p=' num2str(stats.pair(1).signrank.p*(2*lenb))])
                if stats.pair(1).signrank.p < 0.05/(2*lenb)
                    hold on;
%                     p = plot([start(1), (freq(b)+freq(b+1))/2], [1 1]*yl(2), '-', 'color', 0.5*[0 1 0], 'linewidth', 2);
                    p = plot(bandrange{b}, [1 1]*yl(2), '-', 'color', 0.5*[0 1 0], 'linewidth', 2);
                    p.Color(4) = 0.4;
                end
%                 start(1) = (freq(b)+freq(b+1))/2;

                % condition 2
                if ndrug > 1
                    subplot(3*ndrug, 2, 1+2*(c-1)+ndrug*2*(a-1))
                    yl = get(gca, 'YLim');
                    disp(['cond 2:' bandnames{b} ', ' num2str(a) ', p=' num2str(stats.pair(2).signrank.p*(2*lenb))])
                    if stats.pair(2).signrank.p < 0.05/(2*lenb)
                        hold on;
                        p = plot(bandrange{b}, [1 1]*yl(2), '-', 'color', 0.5*[0 1 0], 'linewidth', 2);
                        p.Color(4) = 0.4;
                    end
%                     start(2) = (freq(b)+freq(b+1))/2;
                
%                     % condition 1 vs 2
%                     subplot(3*ndrug, 2, 2+2*(c-1)+ndrug*2*(a-1))
%                     if stats.ranksum.p < 0.05/length(step)-1
%                         hold on;
%                         p = plot(bandrange{b}, [1 1]*(yl(2) - 0.05*(yl(2)-yl(1))), ':', 'color', lcols(d,:), 'linewidth', 2);
%                         p.Color(4) = 0.4;
%                     end
                end
            end
%             for b = 1:lenb
%                 % pair tests
%                 stats = pair_tests(pairs{b, 1}, pairs{b, 2});
% 
%                 % condition 1
%                 subplot(3*ndrug, 2, 2+ndrug*2*(a-1))
%                 if stats.pair(1).signrank.p < 0.05/lenb
%                     hold on;
%                     p = plot(bandrange{b}, [1 1]*yl(2), '-', 'color', lcols(d,:), 'linewidth', 2);
%                     p.Color(4) = 0.4;
%                 end
% 
%                 % condition 2
%                 if ndrug > 1
%                     subplot(3*ndrug, 2, 2+2*(c-1)+ndrug*2*(a-1))
%                     if stats.pair(2).signrank.p < 0.05/lenb
%                         hold on;
%                         p = plot(bandrange{b}, [1 1]*yl(2), '-', 'color', lcols(d,:), 'linewidth', 2);
%                         p.Color(4) = 0.4;
%                     end
% 
%                     % condition 1 vs 2
%                     subplot(3*ndrug, 2, 2+2*(c-1)+ndrug*2*(a-1))
%                     if stats.ranksum.p < 0.05/lenb
%                         hold on;
%                         p = plot(bandrange{b}, [1 1]*(yl(2) - 0.05*(yl(2)-yl(1))), ':', 'color', lcols(d,:), 'linewidth', 2);
%                         p.Color(4) = 0.4;
%                     end
%                 end
%             end
            xlim(xrange)
            ylim(yl)                
        end
    end
    fign = fign + 1;    
    set(gcf, 'Name', fields{f}, 'NumberTitle', 'off')
    savefig(gcf, [mypath ...
        '/Katsuhisa/serotonin_project/LFP_project/Manuscript/Figures/Fig3_IRASA/raw_figs/' ...
        fields{f} '_' splittype '.fig'])
end   

% power-law fits
figure(fign);
for a = 1:3 % animal
    pairs = cell(4, 2);
    for c = 1:ndrug % is5ht
        % condition
        if a < 3
            cond = is5ht==c-1 & ismango==a-1;
        else
            cond = is5ht==c-1;
        end

        for d = 1:2 % base or drug
            % mean power-law fit params
            pairs{1, c}(:, d) = squeeze(para.cond(d).enc(cond, 3, 1));
            pairs{2, c}(:, d) = squeeze(para.cond(d).enc(cond, 4, 1));
            pairs{3, c}(:, d) = squeeze(para.cond(d).enc(cond, 5, 1));
            pairs{4, c}(:, d) = squeeze(para.cond(d).enc(cond, 6, 1));
%             beta0 = squeeze(nanmean(para.cond(d).enc(cond, 3, 1), 1));
%             cons0 = squeeze(nanmean(para.cond(d).enc(cond, 4, 1), 1));
%             beta1 = squeeze(nanmean(para.cond(d).enc(cond, 5, 1), 1));
%             cons1 = squeeze(nanmean(para.cond(d).enc(cond, 6, 1), 1));
            
%             powlaw0 = 10.^(polyval([-beta0 cons0], [freq(1) log10(freq(idx14))]));
%             powlaw1 = 10.^(polyval([-beta1 cons1], log10([freq(idx15) freq(end)])));
%             loglog([freq(1) 14 15 100], [powlaw0, powlaw1], lspc{d}, 'color', lcols(s,:));
        end
        
        % plot 
        subplot(3*ndrug, 2, 1 + 2*(c-1) + 2*ndrug*(a-1))
        ls_scatter(pairs{1, c}(:), pairs{2, c}(:), [zeros(sum(cond==1), 1); ones(sum(cond==1), 1)])
        if a==3 && c==2
            xlabel('power-low exponent (3 - 10 Hz)')
        end
        ylabel({animal_lab{a}, 'intersept', pairnames{c, 2}})
        subplot(3*ndrug, 2, 2 + 2*(c-1) + 2*ndrug*(a-1))
        ls_scatter(pairs{3, c}(:), pairs{4, c}(:), [zeros(sum(cond==1), 1); ones(sum(cond==1), 1)])
        if a==3 && c==2
            xlabel('power-low exponent (15 - 90 Hz)')
        end
        ylabel({animal_lab{a}, 'intersept', pairnames{c, 2}})
    end
    
%     % stats
%     for b = 1:4
%         % pair tests
%         stats = pair_tests(pairs{b, 1}, pairs{b, 2});
% 
%         % condition 1
%         subplot(3*ndrug, 1, 1+ndrug*(a-1))
%         if stats.pair(1).signrank.p < 0.05/lenb
%             disp([animal_lab{a} '; cond 1: ' paranames{1+b} ' significant; p=' num2str(stats.pair(1).signrank.p)])
%         end
% 
%         % condition 2
%         if ndrug > 1
%             subplot(3*ndrug, 1, 2+ndrug*(a-1))
%             if stats.pair(2).signrank.p < 0.05/lenb
%                 disp([animal_lab{a} '; cond 2: ' paranames{1+b} ' significant; p=' num2str(stats.pair(2).signrank.p)])
%             end
% 
%             % condition 1 vs 2
%             subplot(3*ndrug, 1, 2+ndrug*(a-1))
%             if stats.ranksum.p < 0.05/lenb
%                 disp([animal_lab{a} '; cond 1 vs 2: ' paranames{1+b} ' significant; p=' num2str(stats.ranksum.p)])
%             end
%         end
%     end
    fign = fign + 1;
    set(gcf, 'Name', 'power-law fits', 'NumberTitle', 'off')
end  
        
% encoding analysis
if ~strcmp(stmtype, 'rc')
    for k = 1:ndrug % is5ht     
        for j = 1:lenp % activities
            figure(fign+j);
            for e = 1:lene % tuning index
                % plot data
                subplot(ndrug, lene, e + lene*(k - 1))
                try
                    % plot
                    anim = 1*ismango==0;
                    unity_scatter(squeeze(para.cond(1).enc(is5ht==k-1, j, e)), ...
                        squeeze(para.cond(2).enc(is5ht==k-1, j, e)), anim(is5ht==k-1))

                    % format
                    if k==1
                        title(encnames{e})
                        if e==1
                            ylabel({paranames{j}, pairnames{k, 2}})
                        end
                    else
                        xlabel({paranames{j}, pairnames{k, 1}})
                        if e==1
                            ylabel({paranames{j}, pairnames{k, 2}})
                        end
                    end
                catch
                    continue
                end
            end
            set(gcf, 'Name', [animal_lab{a} ': ' paranames{j}], 'NumberTitle', 'off')
        end
    end  
else
    figure(fign+1);
    for k = 1:ndrug % is5ht     
        for j = 1:lenp % activities
            % plot data
            subplot(ndrug, lenp, j + lenp*(k - 1))
            try
                % plot
                anim = 1*ismango==0;
                unity_scatter(squeeze(para.cond(1).enc(is5ht==k-1, j, 1)), ...
                    squeeze(para.cond(2).enc(is5ht==k-1, j, 1)), anim(is5ht==k-1))

                % format
                if k==1
                    title(paranames{j})
                    if j==1
                        ylabel(pairnames{k, 2})
                    end
                else
                    xlabel(pairnames{k, 1})
                    if j==1
                        ylabel(pairnames{k, 2})
                    end
                end
            catch
                continue
            end
            set(gcf, 'Name', [animal_lab{a} ': IRASA params'], 'NumberTitle', 'off')
        end
    end  
end
