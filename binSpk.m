function binSpk(exinfo)
% convert spiking activity into spike count of 10ms bins
% INPUT: exinfo
% OUTPUT: sc ... spike count

savepath = 'Z:\Katsuhisa\LFP_project\Data\SpikeCounts\';
for i = 1:length(exinfo)
    disp(['working on session ' num2str(i)])
    
    % control
    load(exinfo(i).fname)
    spikecount = ex2sc(ex);
    slash = strfind(exinfo(i).fname, '\');
    fname = exinfo(i).fname(slash(end)+1:end);
    save([savepath fname], 'spikecount')
    
    % drug
    load(exinfo(i).fname_drug)
    spikecount = ex2sc(ex);
    slash = strfind(exinfo(i).fname_drug, '\');
    fname = exinfo(i).fname(slash(end)+1:end);
    save([savepath fname], 'spikecount')
end






