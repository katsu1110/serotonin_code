function plot_c2s
%%
% plot results from spike prediction from LFPs
%

if ispc
    mypath = 'Z:/';
else
    mypath = '/gpfs01/nienborg/group/';
end

load([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met.mat'])

metrics = {'correlations', 'info'};
metricnames = {'correlation', 'MI'};
lenm = length(metrics);
ani = {'s', 'o'};
lena = length(ani);
col = {'k', 'r'};
drugs = {'NaCl', '5HT'};

close all;
figure;
for k = 1:lenm
    for d = 1:2
        subplot(lenm, 2, d+2*(k-1))
        x = met(met(:, 2)==d-1, 3+4*(k-1));
        y = met(met(:, 2)==d-1, 5+4*(k-1));
        z = met(met(:, 2)==d-1, 1);
        unity_scatter(x, y, z)
        hold on;
        if k==1
            title([drugs{d}])
        end
        if d==1
            ylabel({metricnames{k}, 'drug'})
        else
            xlabel('baseline')
        end
    end
end