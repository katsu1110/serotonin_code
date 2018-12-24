function interaction_lfp_dataset(dat)
% run 'pupil_interaction.m' and 'pupil_lfp.m' in the RC dataset
% INPUT: load([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat'])

% load data
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group/';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/katsuhisa_analysis']))
addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))
pss = struct('session', []);

% loop for pairs of sessions
lend = length(dat);
out_su = cell(1, lend);
out_mu = cell(1, lend);
% out_lfp = cell(1, lend);
parfor i = 1:lend
    try
        out_su{i} = pupil_interaction(dat(i), 'SU');
    catch
        out_su{i} = nan;        
    end
    try
        out_mu{i} = pupil_interaction(dat(i), 'MU');
    catch
        out_mu{i} = nan;
    end
    try 
        out_lfp{i} = pupil_lfp(dat(i));
    catch
        out_lfp{i} = 0;
    end
%     disp(['------------------ session ' num2str(i) ' is processed! ----------------------'])
end
for i = 1:lend
    if isstruct(out_su{i})
        pss.session(i).psintr_su = out_su{i};
        pss.session(i).psintr_su_exist = 1;
         disp(['done pupil interaction analysis SU at session ' num2str(i)])
    else
        pss.session(i).psintr_su_exist = 0;
        disp(['error in pupil interaction analysis SU at session ' num2str(i)])
    end
    if isstruct(out_mu{i})
        pss.session(i).psintr_mu = out_mu{i};
        pss.session(i).psintr_mu_exist = 1;
        disp(['done pupil interaction analysis MU at session ' num2str(i)])
    else
        pss.session(i).psintr_mu_exist = 0;        
        disp(['error in pupil interaction analysis MU at session ' num2str(i)])
    end
    if isstruct(out_lfp{i})
        pss.session(i).pslfp = out_lfp{i};
        pss.session(i).pslfp_exist = 1;
    else
        pss.session(i).pslfp_exist = 0;        
        disp(['error in pupil vs lfp analysis at session ' num2str(i)])
    end
end
% save data structure
save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/pss.mat'], 'pss','-v7.3') 
disp('output structure pss saved')
