function extract_c2s
%%
% extract info from c2s analysis
%

if ispc
    mypath = 'Z:/';
else
    mypath = '/gpfs01/nienborg/group/';
end

% % load results
% M = csvread([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/c2s/results.csv'], 0, 1);

datapath = [mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/'];
dirs = dir(datapath);
dirs(1:2) = [];
lend = length(dirs);
met = nan(lend, 10); % animal, drug, corr mean (base), corr sem (base), corr mean (drug), corr sem (drug), ...
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
    c = 3;
    for k = 1:length(metrics)
%         try
            data = load([dirs(i).folder '/' dirs(i).name '/data.mat']);
            data = data.data;
            ntr = length(data);
            bd = ones(1, ntr);
            for n = 1:ntr
                bd(n) = data{n}.cell_num;
            end
            out = load([dirs(i).folder '/' dirs(i).name '/' metricnames{k} '.mat']);
            out = out.(metrics{k});
            for b = 1:2
                met(i, c) = nanmean(out(1, bd==b), 2);
                met(i, c + 1) = nanstd(out(1, bd==b), [], 2)/sqrt(sum(bd==b));
                c = c + 2;
            end
%         catch
%             continue
%         end
    end        
end

save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met.mat'], 'met')