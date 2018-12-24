function [serdata, tudata] = addLFP(exinfo, lists, pairtype, internal)
%% run 'LFPanalyzer.m' in batch using exinfo
%
% load('Z:/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat')
% load([mypath 'Corinna/SharedCode/Katsu/list_RC.mat']) <--- lists{1}
% load([mypath 'Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat']) <--- lists{2}
% written by Katsuhisa (29.09.17)
% ++++++++++++++++++++++++++++++++

if nargin < 3; pairtype = 'drug'; end
if nargin < 4; internal = 0; end
serdata = struct('session', []); tudata = struct('session', []);

%%
% stimulus condition
list1 = lists{1};
list2 = lists{2};
lene = length(exinfo);
goodunit = zeros(1, lene);
out = cell(1, lene); out2 = cell(1, lene);
parfor i = 1:lene
    % good unit or not
    if (exinfo(i).isRC==1 && ismember(i, list1)) || ...
         (exinfo(i).isRC==0 && ismember(i, list2))
        goodunit(i) = 1;
    end
    % LFP analysis if the session's stimulus is either 'co' or 'RC'
%     if strcmp(exinfo(i).param1, 'co') || exinfo(i).isRC==1
        try
            [out{i}, out2{i}] = LFPanalyzer(exinfo(i), pairtype, internal);    
            disp(['session ' num2str(i) ': analyzed (pairtype = ' pairtype ')'])
        catch
            out{i} = nan; out2{i} = nan;
            disp(['session ' num2str(i) ': error in LFP analyzer (pairtype = ' pairtype ')'])
        end
%     end
end
for i = 1:lene
    serdata.session(i).goodunit = goodunit(i);
    serdata.session(i).exist = 1;
    serdata.session(i).results = out{i};
    tudata.session(i).goodunit = goodunit(i);
    tudata.session(i).exist = 1;
    tudata.session(i).results = out2{i};
end

%%
% autosave
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group/';
else
    mypath = 'Z:/';
end
save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/serdata_' pairtype num2str(internal)  '.mat'], 'serdata', '-v7.3')
disp(['serdata (' pairtype num2str(internal)  ') was saved!'])
save([mypath 'Katsuhisa/serotonin_project/LFP_project/Data/tudata_' pairtype num2str(internal)  '.mat'], 'tudata', '-v7.3')
disp(['tudata (' pairtype num2str(internal)  ') was saved!'])