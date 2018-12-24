 function anaT = analysis_table(lfps, splittype)
%%
% generate a table to see how variables are related to one another
% sessions x variables x stimulus type
%

if nargin < 2; splittype = 'drug'; end

% single pair or more? ===========================
ndrug = 1;
switch splittype
    case 'drug'
        ndrug = 2;
        pairnames = {'base', 'NaCl'; 'base', '5HT'};
    case 'sc'
        pairnames = {'small sc', 'large sc'; 'small sc', 'large sc'};
    case 'pupil'
        pairnames = {'small ps', 'large ps'; 'small ps', 'large ps'};
end
if isfield(lfps, 'lfplist')
    datast = lfps.LFP_prepro;
    lists = [lfps.goodunit', lfps.is5ht', lfps.animal'];
else
    datast{1} = lfps;
    ndrug = 1;
    lists = [1, lfps.is5ht, lfps.ismango];
end

close all;
disp(['5HT: ' num2str(sum(lists(:,2)==1)) ' pairs'])
disp(['NaCl: ' num2str(sum(lists(:,2)==0)) ' pairs'])
% if ndrug==1 
%     lists(:,2) = 0;
% end
lenses = size(lists, 1);

% exclude no data sessions
if strcmp(datast{1}.stm.param, 'rc')
    outs = zeros(lenses, 1);
%     outs([6:9]) = 1; % mango, unkown noise, c2
%     outs = lists(:,1);
%     outs = 1*(outs < 1);
    for i = 1:lenses
        if datast{i}.cond(1).sta.nspk==0 || datast{i}.cond(2).sta.nspk==0
            outs(i) = 1;
        end
    end
    datast(outs==1) = [];
    lists(outs==1, :) = [];
    lenses = size(lists, 1);
end

% path =======================================
if ismember(1, contains(cd, 'gpfs0'))
    mypath = '/gpfs01/nienborg/group';
elseif ismember(1, contains(cd, '\\172.25.250.112'))
    mypath = '//172.25.250.112/nienborg_group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/CircStat2012a']))
addpath(genpath([mypath '/Katsuhisa/serotonin_project/chronux_2_12/chronux_2_12']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/flexibleModulatedPoisson/code_flexibleModulatedPoisson']))

% load structures
load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/smlinfo.mat'])
load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/SI/SI_all.mat'])

% stimulus type ============================
stmdur = 0.45;
dur = 0.2;
switch datast{1}.stm.param
    case 'rc'
%         % remove outliers
%         datast(1:3) = [];
%         lists(1:3, :) = [];
        
        stmidx = ones(lenses, 1);
        stmlab = {'ORxRC'};
        stmdur = 2;
        dur = 1.2;
    case 'or'
        stmidx = zeros(lenses, 2);
        stmlab = {'anti-pref', 'pref'};
        for i = 1:lenses
            [~, sortidx] = sort(datast{i}.cond(1).spk_tu{2}(:,1));
            mini = 2;
            maxi = length(sortidx);
            try
                while isnan(datast{i}.cond(1).sta.p{sortidx(mini)})...
                        | isnan(datast{i}.cond(2).sta.p{sortidx(mini)})
                    mini = mini + 1;
                end
                while isnan(datast{i}.cond(1).sta.p{sortidx(maxi)})...
                        | isnan(datast{i}.cond(2).sta.p{sortidx(maxi)})
                    maxi = maxi - 1;
                end
                stmidx(i, :) = [sortidx(mini), sortidx(maxi)];
            catch
                disp('')
            end
        end
    case 'co'
        stmidx = zeros(lenses, 3);
        stmlab = {'0.25', '0.5', '1'};
        cos = [0.25, 0.5, 1];
        for i = 1:lenses
            for c = 1:3
                [~, idx] = min(abs(datast{i}.stm.vals - cos(c)));
                stmidx(i, c) = idx;
            end
        end
    case {'sz', 'sf'} 
        stmidx = zeros(lenses, 2);
        stmlab = {'min res stm', 'max res stm'};
        for i = 1:lenses
            [~, sortidx] = sort(datast{i}.cond(1).spk_tu{2}(:,1));
            mini = sortidx(1);
            if mini == 1
                mini = sortidx(2);
            end
            stmidx(i, :) = [mini, sortidx(end)];
        end
end

% for oher structs
lfplist = lfps.lfplist(outs==0);
lenses = size(lists, 1);
sesidx = nan(lenses, 1);
lensm = length(smlinfo.fname);
for i = 1:lenses
    for k = 1:lensm
        if contains(smlinfo.fname{k}, lfplist{i}{1})
            sesidx(i) = k;
        end
    end
end
datast = encoding_tuning4rc(datast);

% remove no-data sessions
datast(stmidx(:,1)==0) = [];
lists(stmidx(:,1)==0, :) = [];
stmidx(stmidx(:,1)==0, :) = [];

% preset params
ss = size(stmidx, 2);

% variable names
varnames = {'fr base', 'fr drug', 'r rate', 'noise corr base', 'noise corr drug', 'd noise corr', ...
    'selectivity base', 'selectivity drug', 'd selectivity', 'snr2 base', 'snr2 drug', 'd snr2', ...
    'tuning err base', 'tuning err drug', 'd tuning err', 'sig corr base', 'sig corr drug', 'd sig corr', 'ff base', 'ff drug', 'd ff',...
    'median si-mu base', 'median si-mu drug', 'd median si-mu', 'var si-mu base', 'var si-mu drug', 'd var si-mu', ...
    'additive change', 'gain change', 'wavewidth', 'RFx', 'RFy', 'sta amp base', 'sta amp drug', ...
    'd sta amp', 'sta Tmin base', 'sta Tmin drug', 'd sta Tmin', 'theta pow base', 'theta pow drug', ...
    'd theta pow', 'alpha pow base', 'alpha pow drug', 'd alpha pow', 'beta pow base', 'beta pow drug', ...
    'd beta pow', 'gamma pow base', 'gamma pow drug', 'd gamma pow', 'low-freq pow base', 'low-freq pow drug', 'd low-freq pow', ...
    'broadband pow base', 'broadband pow drug', 'd broadband pow', 'theta res base', 'theta res drug', 'd theta res', ...
    'alpha res base', 'alpha res drug', 'd alpha res', 'beta res base', 'beta res drug', ...
    'd beta res', 'gamma res base', 'gamma res drug', 'd gamma res', 'low-freq res base', 'low-freq res drug', 'd low-freq res', ...
    'broadband res base', 'broadband res drug', 'd broadband res', 'coh theta base', 'coh theta drug', 'd coh theta', ...
    'coh alpha base', 'coh alpha drug', 'd coh alpha', 'coh beta base', 'coh beta drug', 'd coh beta', ...
    'coh gamma base', 'coh gamma drug', 'd coh gamma','coh low-freq base', 'coh low-freq drug', 'd coh low-freq',...
    'coh broadband base', 'coh broadband drug', 'd coh broadband', 'phase theta base', 'phase theta drug', 'd phase theta', ...
    'phase alpha base', 'phase alpha drug', 'd phase alpha', 'phase beta base', 'phase beta drug', 'd phase beta', ...
    'phase gamma base', 'phase gamma drug', 'd phase gamma','phase low-freq base', 'phase low-freq drug', 'd phase low-freq',...
    'phase broadband base', 'phase broadband drug', 'd phase broadband', 'lfpsd theta base', 'lfpsd theta drug', 'd lfpsd theta', ...
    'lfpsd alpha base', 'lfpsd alpha drug', 'd lfpsd alpha', 'lfpsd beta base', 'lfpsd beta drug', 'd lfpsd beta', ...
    'lfpsd gamma base', 'lfpsd gamma drug', 'd lfpsd gamma', 'lfpsd low-freq base', 'lfpsd low-freq drug', ...
    'd lfpsd low-freq', 'lfpsd broadband base', 'lfpsd broadband drug', 'd lfpsd broadband', ...
    'lfpsd signal base', 'lfpsd signal drug', 'd lfpsd signal', 'lfppow ratio base', 'lfppow ratio drug', 'd lfppow ratio'};

lenv = length(varnames);

bandnames = {'theta (3-7)', 'alpha (8-13)', 'beta (14-24)', 'gamma (36-48)', ...
    'low-freq (3-10)', 'broadband (3-48)'};
bandrange = {[3, 7], [8, 12], [15, 25], [40, 48], [3, 10], [3, 48]};
lenb = length(bandnames);
    
at = nan(lenses, lenv, ss);

% bandpass filtering for spike counts ===============================
% [b1, a1] = butter(8, 1/200, 'high');

for i = 1:lenses
    for s = 1:ss
        % firing rate and sensory coding =====================
        % fr base
        at(i, 1, s) = datast{i}.cond(1).spk_tu{1}(1);
        % fr drug
        at(i, 2, s) = datast{i}.cond(2).spk_tu{1}(1);
        % r rate
        at(i, 3, s) = datast{i}.cond(2).spk_tu{1}(s, 1)/datast{i}.cond(1).spk_tu{1}(s, 1);
        % noise corr base
        sua0ori = round(dur*datast{i}.cond(1).lfprel.mat{s}(:, 1))';
%         sua0 = filter(b1, a1, sua0ori);
        sua0 = locdetrend(sua0ori, 1, [20, 1]);
%         sua0 = sua0 - mean(sua0);
        sua0 = sua0 - mean(sua0) + dur*at(i, 1, s);
        sua0(sua0 <= 0) = 0;
        sua0 = round(sua0);
        mua0 = round(dur*datast{i}.cond(1).lfprel.mat{s}(:, 2))';
%         mua0 = filter(b1, a1, mua0);
        mua0 = locdetrend(mua0, 1, [20, 1]);
%         mua0 = mua0 - mean(mua0);
        mua0 = round(mua0 - mean(mua0) + dur*mean(datast{i}.cond(1).lfprel.mat{s}(:, 2)));
        rr = corrcoef(zscore(sua0), zscore(mua0));
        at(i, 4, s) = rr(1,2);
        at(i, 4, s) = 0.5*log((1+at(i, 4, s))/(1-at(i, 4, s)));
        % noise corr drug
        sua2ori = dur*datast{i}.cond(2).lfprel.mat{s}(:, 1)';
%         sua2 = filter(b1, a1, sua2ori);
        sua2 = locdetrend(sua2ori, 1, [20, 1]);
%         sua2 = sua2 - mean(sua2);
        sua2 = sua2 - mean(sua2) + dur*at(i, 2, s);
        sua2(sua2 <= 0) = 0;
        sua2 = round(sua2);
        mua2 = dur*datast{i}.cond(2).lfprel.mat{s}(:, 2)';
%         mua2 = filter(b1, a1, mua2);
        mua2 = locdetrend(mua2, 1, [20, 1]);
%         mua2 = mua2 - mean(mua2);
        mua2 = round(mua2 - mean(mua2) + dur*mean(datast{i}.cond(2).lfprel.mat{s}(:, 2)));
        rr = corrcoef(zscore(sua2), zscore(mua2));
        at(i, 5, s) = rr(1,2);
        at(i, 5, s) = 0.5*log((1+at(i, 5, s))/(1-at(i, 5, s)));
        % d noise corr
        at(i, 6, s) = at(i, 4, s) - at(i, 5, s);
        % selectivity base
        at(i, 7, s) = datast{i}.cond(1).tuning.spikes.unique.circularvariance;
        % selectivity drug
        at(i, 8, s) = datast{i}.cond(2).tuning.spikes.unique.circularvariance;
        % d selectivity
        at(i, 9, s) = at(i, 7, s) - at(i, 8, s);
        % snr2 base
        at(i, 10, s) = nanmean(datast{i}.cond(1).tuning.spikes.snr2);
        % snr2 drug
        at(i, 11, s) = nanmean(datast{i}.cond(2).tuning.spikes.snr2);
        % d snr2
        at(i, 12, s) = at(i, 10, s) - at(i, 11, s);
        % tuning error base
        at(i, 13, s) = nanmean(datast{i}.cond(1).tuning.su_fitparam.val.sem);
        % tuning error drug
        at(i, 14, s) = nanmean(datast{i}.cond(2).tuning.su_fitparam.val.sem);
        % tuning error base
        at(i, 15, s) = at(i, 13, s) - at(i, 14, s);
        % signal corr base
        rr = corrcoef(zscore(datast{i}.cond(1).tuning.su_fitparam.val.mn), ...
            zscore(datast{i}.cond(1).tuning.mu_fitparam.val.mn));
        at(i, 16, s) = rr(1, 2);
        % signal corr drug
        rr = corrcoef(zscore(datast{i}.cond(2).tuning.su_fitparam.val.mn), ...
            zscore(datast{i}.cond(2).tuning.mu_fitparam.val.mn));
        at(i, 17, s) = rr(1, 2);
        % d signal corr
        at(i, 18, s) = at(i, 16, s) - at(i, 17, s);
        % fano factor base
        prshat_srp = soft_rect_p(sua0);
        at(i, 19, s) = exp(prshat_srp(end))/prshat_srp(1);
%          at(i, 19, s) = std(sua0);
%         at(i, 19, s) = var(sua0)/mean(sua0);
%         at(i, 19, s) = sqrt(exp(var(log(sua0))) - 1);
        % fano factor drug
        prshat_srp = soft_rect_p(sua2);
        at(i, 20, s) = exp(prshat_srp(end))/prshat_srp(1);
%         at(i, 20, s) = var(sua2)/mean(sua2);
%         at(i, 20, s) = std(sua2);
%         at(i, 20, s) = sqrt(exp(var(log(sua2))) - 1);
        % d fano factor
        at(i, 21, s) = at(i, 19, s) - at(i, 20, s);
        
        % within-trial variability ===========================
        % median si base
        at(i, 22, s) = nanmedian(SI_all.si{1}{sesidx(i)}(:, 2));
        % median si drug
        at(i, 23, s) = nanmedian(SI_all.si{2}{sesidx(i)}(:, 2));
        % d median si
        at(i, 24, s) = at(i, 22, s) - at(i, 23, s);
        % var si base
        at(i, 25, s) = nanvar(SI_all.si{1}{sesidx(i)}(:, 2));
        % var si drug
        at(i, 26, s) = nanvar(SI_all.si{2}{sesidx(i)}(:, 2));
        % d var si
        at(i, 27, s) = at(i, 25, s) - at(i, 26, s);
        
        % change in tuning curve ====================
        beta = gmregress(datast{i}.cond(1).tuning.su_fitparam.val.mn, ...
            datast{i}.cond(2).tuning.su_fitparam.val.mn);
        % additive change
        at(i, 28, s) = beta(1);
        % gain change 
        at(i, 29, s) = beta(2);
        
        % neural properties ======================
        % wavewidth
        at(i, 30, s) = smlinfo.paramat(sesidx(i), 7);
        % RF x
        at(i, 31, s) = smlinfo.paramat(sesidx(i), 8);
        % RF y
        at(i, 32, s) = smlinfo.paramat(sesidx(i), 9);
        
        % stLFP ================================
        % sta amp base
        sta_t = linspace(-datast{i}.wnd, datast{i}.wnd, ...
            size(datast{i}.cond(1).sta.mean, 2));
%         sta_analysiswnd = sta_t >0 & sta_t < 0.03;
        sta_analysiswnd = sta_t >= -0.05 & sta_t <= 0.05;
        sta_analysis_time = sta_t(sta_analysiswnd);
        [maxv, maxt] = max(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));
        [minv, zerotb] = min(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_analysiswnd));
%         [minv, zerotb] = min(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_t < maxt));
%         reft = sta_t(sta_t < maxt);
%         zerotb = reft(zerotb);
        zerotb = sta_analysis_time(zerotb);
%         at(i, 30, s) = abs(ampb);
        at(i, 33, s) = maxv - minv;
        % sta amp drug
        [maxv, maxt] = max(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));
        [minv, zerotd] = min(datast{i}.cond(2).sta.mean(stmidx(i, s), sta_analysiswnd));
%         [minv, zerotb] = min(datast{i}.cond(1).sta.mean(stmidx(i, s), sta_t < maxt));
%         reft = sta_t(sta_t < maxt);
%         zerotb = reft(zerotb);
        zerotd = sta_analysis_time(zerotd);
%         at(i, 31, s) = abs(ampd);
        at(i, 34, s) = maxv - minv;
        % d sta amp
        at(i, 35, s) = at(i, 33, s) - at(i, 34, s);
        % sta Tmin base
        at(i, 36, s) = zerotb;
        % sta Tmin drug
        at(i, 37, s) = zerotd;
        % d sta Tmin
        at(i, 38, s) = at(i, 36, s) - at(i, 37, s);
        
        % spectral power & response & coherence & phase & lfp SD & LFP vs PS or FR ======
        freq = datast{i}.cond(1).spectrogram.f{1};
        spe_t = datast{i}.cond(1).spectrogram.t{1};
        lent = length(spe_t);
        spe_t = linspace(-0.1, stmdur, lent);
        S1 = 10*log10(datast{i}.cond(1).spectrogram.S{stmidx(i, s)});
        S2 = 10*log10(datast{i}.cond(2).spectrogram.S{stmidx(i, s)});
        freqc = datast{i}.cond(1).coherence.f{stmidx(i, s)};
        coh1 = datast{i}.cond(1).coherence.C{stmidx(i, s)};
        coh2 = datast{i}.cond(2).coherence.C{stmidx(i, s)};
        phase1 = datast{i}.cond(1).coherence.phi{stmidx(i, s)};
        phase2 = datast{i}.cond(2).coherence.phi{stmidx(i, s)};
        for b = 1:lenb
            % frequency range
            frange = bandrange{b}(1) <= freq & bandrange{b}(2) >= freq;
            frangec = bandrange{b}(1) <= freqc & bandrange{b}(2) >= freqc;
            
            % overall power =============
            % base 
            at(i, 39 + 3*(b-1), s) = nanmean(nanmean(S1(frange, spe_t > datast{i}.window{end}(1))));
            % drug 
            at(i, 40 + 3*(b-1), s) = nanmean(nanmean(S2(frange, spe_t > datast{i}.window{end}(1))));
            % d
            at(i, 41 + 3*(b-1), s) = at(i, 39 + 3*(b-1), s) - at(i, 40 + 3*(b-1), s);
            
            % response power ===========
            % base 
            at(i, 57 + 3*(b-1), s) = at(i, 39 + 3*(b-1), s) - ...
                nanmean(nanmean(S1(frange, spe_t < 0)));
            % drug 
            at(i, 58 + 3*(b-1), s) = at(i, 40 + 3*(b-1), s) - ...
                nanmean(nanmean(S2(frange, spe_t < 0)));
            % d
            at(i, 59 + 3*(b-1), s) = at(i, 57 + 3*(b-1), s) - at(i, 58 + 3*(b-1), s);            
            
            % coherence ===============
            % base 
            at(i, 75 + 3*(b-1), s) = nanmean(coh1(frangec));
            % drug
            at(i, 76 + 3*(b-1), s) = nanmean(coh2(frangec));
            % d
            at(i, 77 + 3*(b-1), s) = at(i, 75 + 3*(b-1), s) - at(i, 76 + 3*(b-1), s);     
            
            % phase ==================
            % base 
            at(i, 93 + 3*(b-1), s) = maxcoh_rad(coh1(frangec), phase1(frangec));
            % drug
            at(i, 94 + 3*(b-1), s) = maxcoh_rad(coh2(frangec), phase2(frangec));
            % d
            at(i, 95 + 3*(b-1), s) = at(i, 93 + 3*(b-1), s) - at(i, 94 + 3*(b-1), s);     
                        
            % LFP SD ==============
            % base 
            at(i, 111 + 3*(b-1), s) = nanstd(datast{i}.cond(1).lfprel.mat{s}(:, 6+b));
            % drug
            at(i, 112 + 3*(b-1), s) = nanstd(datast{i}.cond(2).lfprel.mat{s}(:, 6+b));
            % d
            at(i, 113 + 3*(b-1), s) = at(i, 111 + 3*(b-1), s) - at(i, 112 + 3*(b-1), s);   
        end  
        % LFP SD base 
        at(i, 129, s) = nanstd(datast{i}.cond(1).lfprel.mat{s}(:, 5));
        % LFP SD drug
        at(i, 130, s) = nanstd(datast{i}.cond(2).lfprel.mat{s}(:, 5));
        % d LFP SD
        at(i, 131, s) = at(i, 129, s) - at(i, 130, s);  
        % LFP power ratio base 
        at(i, 132, s) = at(i, 39 + 3*(1-1), s)/at(i, 39 + 3*(4-1), s);
        % LFP power ratio drug
        at(i, 133, s) = at(i, 40 + 3*(1-1), s)/at(i, 40 + 3*(4-1), s);
        % d LFP power ratio
        at(i, 134, s) = at(i, 132, s) - at(i, 133, s);  
    end    
end
% fillnan
for i = 1:size(at, 2)
    nans = isnan(at(:, i));
    at(nans, i) = nanmedian(at(:, i));
end

% autosave
anaT.table = at;
anaT.varnames = varnames;
anaT.lists = lists;
anaT.pairnames = pairnames;
save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/tables/anaT_' splittype '.mat'], 'anaT')
disp('table saved!')

% subfunction
function datast = encoding_tuning4rc(datast)
lenses = length(datast);
for i = 1:lenses
    for d = 1:2
        % tuning values
        or = datast{i}.cond(d).tuning.su_fitparam.val.or;
        mn = datast{i}.cond(d).tuning.su_fitparam.val.mn;
        sem = datast{i}.cond(d).tuning.su_fitparam.val.sem;
        
        % mean
        datast{i}.cond(d).tuning.spikes.mean = mn;
        
        % SNR2
        snr2 = (mn./sem).^2;
        snr2(isnan(snr2)) = 0;
        snr2(isinf(snr2)) = 100;
        datast{i}.cond(d).tuning.spikes.snr2 = snr2;
        
        % circular variance --- Ringach et al. (2002)
        or = or * pi /180;
        or = mod(or, 2*pi);
        if sum(mn < 0) > 0
            mn = mn + abs(min(mn));
        end
        % compute weighted sum of cos and sin of angles
        r = sum(mn.*exp(1i*or));
        % obtain length 
        r = abs(r)./sum(mn);
        datast{i}.cond(d).tuning.spikes.unique.circularvariance = 1 - r;
    end
end

function rad = maxcoh_rad(coh, rad)
lenv = length(coh);
if lenv < 4
    rad = circ_mean(rad);
else
    [~, maxi] = max(coh);
    if maxi==1
        idx = [1 2 3];
    elseif maxi==lenv
        idx = [lenv-2 lenv-1 lenv];
    else
        idx = [maxi-1 maxi maxi+1];
    end
    rad = circ_mean(rad(idx));
end

