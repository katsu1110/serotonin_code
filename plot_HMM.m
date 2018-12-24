function plot_HMM(hmms)
%%
% visualize results from hmm analysis
%

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))


% data extraction ===================
load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat'], 'dat')
lend = length(dat);
lists = zeros(lend, 2); % is5ht, ismango
for i = 1:lend
    % session info
    lists(i, 1) = dat(i).is5HT;
    if strcmp(dat(i).monkey, 'ma')
        lists(i, 2) = 1;
    end
    
    
end 


% firing rate