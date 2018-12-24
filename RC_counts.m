function RC_counts

load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\LFPprepro\lfplist.mat')
c = zeros(2, 2);
for i = 1:length(lfplist)
    fname = lfplist{i}{2};
    if isnan(fname)
        continue
    end
    fname
    if contains(fname, 'xRC') && contains(fname(1:2), 'ma') && contains(fname, 'NaCl')
        c(1, 1) = c(1, 1) + 1;
    elseif contains(fname, 'xRC') && contains(fname(1:2), 'ka') && contains(fname, 'NaCl')
        c(1, 2) = c(1, 2) + 1;
    elseif contains(fname, 'xRC') && contains(fname(1:2), 'ma') && contains(fname, '5HT')
        c(2, 1) = c(2, 1) + 1;
    elseif contains(fname, 'xRC') && contains(fname(1:2), 'ka') && contains(fname, '5HT')
        c(2, 2) = c(2, 2) + 1;
    end
end
c