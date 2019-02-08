function fix_dim_data

mypath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data';

listings = dir(mypath);
for i = 1:length(listings)
    try
        load([listings(i).folder '/' listings(i).name '/data.mat'])
        for j = 1:length(data)
            data{j}.spike_times = data{j}.spike_times';
        end
        save([listings(i).folder '/' listings(i).name '/data.mat'], 'data')
        disp([listings(i).folder '/' listings(i).name ' fixed and saved!'])
    catch
        continue
    end
end

