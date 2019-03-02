function plot_c2s
%%
% plot results from spike prediction from LFPs
%

if ispc
    mypath = 'Z:/';
else
    mypath = '/gpfs01/nienborg/group/';
end

% load([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met.mat'])
load([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/c2s/met_cv10.mat'])

metrics = {'correlations', 'info'};
metricnames = {'correlation', 'MI'};
lenm = length(metrics);
ani = {'s', 'o'};
lena = length(ani);
col = {'k', 'r'};
drugs = {'NaCl', '5HT'};

close all;
figure;
for k = 1:lenm % metric
    for d = 1:2 % NaCl or 5HT
        subplot(lenm, 3, d+3*(k-1))
        x = met(met(:, 2)==d-1, 3+2*(k-1));
        y = met(met(:, 2)==d-1, 4+2*(k-1));
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
    % FR control
    subplot(lenm, 3, 3+3*(k-1))
    x = met(:, 7+2*(k-1));
    y = met(:, 8+2*(k-1));
    z = met(:, 1);
    unity_scatter(x, y, z)
    hold on;
    if k==1
        title('FR control')
    else
        xlabel('high FR')
    end
    ylabel('low FR')
end