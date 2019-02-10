function format_c2s_eval

mypath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data';

listings = dir(mypath);
fnames = {'train', 'test'};
for i = 1:length(listings)
    for f = 1:length(fnames)
        try
            load([listings(i).folder '/' listings(i).name '/preprocessed_' fnames{f} '.mat'])
            cell_num = ones(1, length(data));
            for j = 1:length(data)
                cell_num(j) = data{j}.cell_num;
            end
            data = data(cell_num==2);
            save([listings(i).folder '/' listings(i).name '/preprocessed_' fnames{f} '.mat'], 'data')
            load([listings(i).folder '/' listings(i).name '/predicted_' fnames{f} '.mat'])
            data = data(cell_num==2);
            save([listings(i).folder '/' listings(i).name '/predicted_' fnames{f} '.mat'], 'data')
            disp([listings(i).folder '/' listings(i).name ' fixed and saved!'])
        catch
            continue
        end
    end
end