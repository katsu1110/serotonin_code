function fix_dim_data

mypath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data';

listings = dir(mypath);
listings = listings(3:end);
for i = 1:length(listings)
%     try
        load([listings(i).folder '/' listings(i).name '/data.mat'])
        cell_num = ones(1, length(data));
        for j = 1:length(data)
%             data{j}.spike_times = data{j}.spike_times';
            cell_num(j) = data{j}.cell_num;
        end
        uni = unique(cell_num);
        
%         % KK's initial analysis
%         % ==========================================================
%         % baseline
%         data_base = data(cell_num==uni(1));
%         [data_base, data1, data2] = split_traintest(data_base);        
%         data_drug = data(cell_num==uni(2));
%         data = data_base;
%         data = data2vec(data);
%         save([listings(i).folder '/' listings(i).name '/data_base.mat'], 'data')
%         
%         % baseline 1 + drug
%         data = data_concatenate(data1, data_drug);
%         data = data2vec(data);
%         save([listings(i).folder '/' listings(i).name '/data_train.mat'], 'data')
%         
%         % baseline 2 + drug
%         data = data_concatenate(data2, data_drug);
%         data = data2vec(data);
%         save([listings(i).folder '/' listings(i).name '/data_test.mat'], 'data')
%         
%         % FR control
%         data_fr = split_traintest(data_base, 'FRcontrol');
%         data = data2vec(data_fr);
%         save([listings(i).folder '/' listings(i).name '/data_fr.mat'], 'data')
%         
%         % =============================================================
        
        % what HN thinks good
        % ==========================================================
        data_ori = data;
        
        % baseline 
        data_base = data_ori(cell_num==uni(1));
        data_base = split_traintest(data_base, 'random', 10);        
        data = data2vec(data_base);
        save([listings(i).folder '/' listings(i).name '/data_base_cv10.mat'], 'data')
        
        % drug
        data_drug = data_ori(cell_num==uni(2));
        data_drug = split_traintest(data_drug, 'random', 10);        
        data = data2vec(data_drug);
        save([listings(i).folder '/' listings(i).name '/data_drug_cv10.mat'], 'data')
                        
        % FR control
        [~, data1, data2] = split_traintest(data_base, 'FRcontrol');
        data1 = split_traintest(data1, 'random', 10);
        data = data2vec(data1);
        save([listings(i).folder '/' listings(i).name '/data_lowFR_cv10.mat'], 'data')
        data2 = split_traintest(data2, 'random', 10);
        data = data2vec(data2);
        save([listings(i).folder '/' listings(i).name '/data_highFR_cv10.mat'], 'data')
        
        % =============================================================
        
        disp([listings(i).folder '/' listings(i).name ' fixed and saved!'])
%     catch
%         disp([listings(i).folder '/' listings(i).name ' ERR!'])
%         continue
%     end
end

function data = data2vec(olddata)
lend = length(olddata);
cell_num = ones(1, lend);
for i = 1:lend
    cell_num(i) = olddata{i}.cell_num;
end
lenc = length(unique(cell_num));
data = cell(1, lenc);
for i = 1:lenc
    data{i}.fps = olddata{1}.fps;
    data{i}.cell_num = i;
    data{i}.spike_times = [];
    data{i}.calcium = [];
end
c = zeros(1, lenc);
for i = 1:lend
    % spike times
    spkt = olddata{i}.spike_times + c(olddata{i}.cell_num);
    data{olddata{i}.cell_num}.spike_times = [data{olddata{i}.cell_num}.spike_times, ...
        spkt];
    
    % LFP
    data{olddata{i}.cell_num}.calcium = [data{olddata{i}.cell_num}.calcium, ...
        olddata{i}.calcium];
    
    c(olddata{i}.cell_num) = c(olddata{i}.cell_num) + 1000*length(olddata{i}.calcium)/data{olddata{i}.cell_num}.fps;
end


function data = data_concatenate(data1, data2)
len1 = length(data1);
len2 = length(data2);
data = cell(1, len1+len2);
c = 1;
for i = 1:len1
    data{c} = data1{i};
    data{c}.cell_num = 1;
    c = c + 1;
end
for i = 1:len2
    data{c} = data2{i};
    data{c}.cell_num = 2;
    c = c + 1;
end

function [data, data1, data2] = split_traintest(data, splittype, n_split)
if nargin < 2; splittype = 'random'; end
if nargin < 3; n_split = 2; end
lend = length(data);
switch splittype
    case 'random'
        rng(19891220);
        idx = randi(n_split, 1, lend);
        for i = 1:lend
            data{i}.cell_num = idx(i);
        end
    case 'FRcontrol'
        spikecounts = zeros(1, lend);
        for i = 1:lend
            spikecounts(i) = length(data{i}.spike_times);
        end
        [~, sortidx] = sort(spikecounts, 'ascend');
        
        % median split
        idx = zeros(1, lend);
        for i = 1:lend
            if i < lend/2
                idx(sortidx(i)) = 1;
                data{sortidx(i)}.cell_num = 1;
            else
                idx(sortidx(i)) = 2;
                data{sortidx(i)}.cell_num = 2;
            end
        end
end
data1 = data;
data1 = data1(idx==1);
data2 = data;
data2 = data2(idx==2);


% function [data_train, data_test] = split_traintest(data, weights)
% rng(19891220);
% lend = length(data);
% idx = datasample([1 2], lend, 'Weights', weights);
% data_train = data(idx==1);
% data_test = data(idx==2);
% idx = randi(2, 1, lend);
% for i = 1:lend
%     data{i}.cell_num = idx(i);
% end