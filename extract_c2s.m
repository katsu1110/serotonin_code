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

% animal, drug, corr (base), corr (drug), MI (base), MI (drug), corr (hFR), corr (lFR), MI (hFR), MI (lFR)
met = nan(lend, 10); 
mfnames = {'correlation', 'MI'};
fieldnames = {'correlations', 'info'};
lenm = length(mfnames);
for i = 1:length(dirs)
    disp(['ses ' num2str(i) ': ' dirs(i).name])
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
%         % baseline
%         data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_base.mat']);
%         met(i, 3+2*(m-1)) = nanmean(data.(fieldnames{m})(2, :), 2);
        
%         % drug
%         data1 = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_train.mat']);
%         data2 = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_test.mat']);
%         data.(fieldnames{m}) = [data1.(fieldnames{m}), data2.(fieldnames{m})];
%         met(i, 4+2*(m-1)) = nanmean(data.(fieldnames{m})(2, :), 2);
        
%         % FR control
%         data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_fr.mat']);
%         met(i, 7+2*(m-1)) = data.(fieldnames{m})(2, 2);
%         met(i, 8+2*(m-1)) = data.(fieldnames{m})(2, 1);
        
        % baseline
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_base_cv10.mat']);
        met(i, 3+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        
        % drug
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_drug_cv10.mat']);
        met(i, 4+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        
        % FR control
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_lowFR_cv10.mat']);
        met(i, 7+2*(m-1)) = model_perf(data, 1, fieldnames{m});
        data = load([dirs(i).folder '/' dirs(i).name '/' mfnames{m} '_highFR_cv10.mat']);
        met(i, 8+2*(m-1)) = model_perf(data, 1, fieldnames{m});
    end
end

% save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met.mat'], 'met')
save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met_cv10.mat'], 'met')

function mpf = model_perf(data, row, name)
if strcmp(name, 'info')
    mpf = nanmean(data.info(row, :)./data.entropy(row, :), 2);
else
    mpf = nanmean(data.(name)(row, :), 2);
end