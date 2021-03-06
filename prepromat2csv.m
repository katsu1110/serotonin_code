function prepromat2csv
%%
% to use deep nets convert 'preprocessed' data by c2s into csv


mypath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data';

listings = dir(mypath);
listings = listings(3:end);
fnames = {'preprocessed_base_cv10', 'preprocessed_drug_cv10', ...
    'preprocessed_highFR_cv10', 'preprocessed_lowFR_cv10'};
delnames = {'auroc.mat', 'correlation.mat', 'correlation_base.mat', ...
    'correlation_fr.mat', 'correlation_test.mat', 'correlation_train.mat', 'cvscores.csv', ...
    'data_drug.pck', 'data_drug.pickle', 'data_test', 'data_test.pck', 'data_test.pickle', ...
    'data_train.pck', 'data_train.pickle', 'MI.mat', 'MI_base.mat', 'MI_fr.mat', 'MI_test.mat', ...
    'MI_train.mat', 'predictions.mat', 'preprocessed.mat', 'preprocessed_base.mat', ...
    'preprocessed_fr.mat', 'preprocessed_test.mat', 'preprocessed_train.mat', 'stlfp0.mat', 'stlfp1.mat', ...
    'data_train.mat', 'data_test.mat', 'predicted_train.mat', 'predicted_test.mat', 'predicted_train2.mat', 'predicted_test2.mat', ...
    'preprocessed2_train.mat', 'preprocessed2_test.mat'};
for i = 1:length(listings)
    % save as csv
    for f = 1:4
        load([listings(i).folder '/' listings(i).name '/' fnames{f} '.mat'])
        lens = zeros(1, 10);
        for k = 1:10
            lens(k) = length(data{k}.spikes);
        end
        max_len = max(lens);
        lfp = nan(max_len, 10); spk = nan(max_len, 10);
        for k = 1:10
            lfp(1:lens(k), k) = data{k}.calcium;
            spk(1:lens(k), k) = data{k}.spikes;
        end
        dlmwrite([listings(i).folder '/' listings(i).name '/lfp' fnames{f}(13:end) '.csv'], lfp, 'delimiter', ',', 'precision', 9)
        dlmwrite([listings(i).folder '/' listings(i).name '/spk' fnames{f}(13:end) '.csv'], spk, 'delimiter', ',', 'precision', 9)
    end
    disp([listings(i).folder '/' listings(i).name ' ... csv saved!'])
    
    % delete unnecessary files
    for d = 1:length(delnames)
        if exist([listings(i).folder '/' listings(i).name '/' delnames{d}], 'file')==2
            delete([listings(i).folder '/' listings(i).name '/' delnames{d}]);
        end
    end
end