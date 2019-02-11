function extract_c2s
%%
% extract info from c2s analysis
%

if ispc
    mypath = 'Z:/';
else
    mypath = '/gpfs01/nienborg/group/';
end

datapath = [mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/'];
dirs = dir(datapath);
dirs(1:2) = [];
lend = length(dirs);
met = nan(lend, 10); % animal, drug, corr mean (base), corr sem (base), corr mean (drug), corr sem (drug), ...
mfnames = {'correlation', 'MI'};
fieldnames = {'correlations', 'info'};
lenm = length(mfnames);
for i = 1:length(dirs)
    if contains(dirs(i).name, 'ka_')
        met(i,1) = 0;
    elseif contains(dirs(i).name, 'ma_')
        met(i, 1) = 1;
    end
    if contains(dirs(i).name, '5HT')
        met(i,2) = 1;
    elseif contains(dirs(i).name, 'NaCl')
        met(i, 2) = 0;
    end
    
    % performance
    for m = 1:lenm
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_base.mat']);
        met(i, 3+4*(m-1)) = nanmean(data.(fieldnames{m})(2, :), 2);
        met(i, 4+4*(m-1)) = nanstd(data.(fieldnames{m})(2, :), [], 2)/sqrt(size(data.(fieldnames{m}), 2));
        data1 = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_train.mat']);
        data2 = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_test.mat']);
        data.(fieldnames{m}) = [data1.(fieldnames{m}), data2.(fieldnames{m})];
        met(i, 5+4*(m-1)) = nanmean(data.(fieldnames{m})(2, :), 2);
        met(i, 6+4*(m-1)) = nanstd(data.(fieldnames{m})(2, :), [], 2)/sqrt(size(data.(fieldnames{m}), 2));
    end
end

save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met.mat'], 'met')