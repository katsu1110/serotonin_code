function plot_c2s
%%
% plot results from spike prediction from LFPs
%

if ispc
    mypath = 'Z:';
else
    mypath = '/gpfs01/nienborg/group';
end

% load results
M = csvread([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/c2s/results.csv'], 0, 1);