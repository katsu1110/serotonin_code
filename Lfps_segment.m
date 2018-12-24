function Lfps = Lfps_segment(Lfps, sesnum)
% smaller lfps

fieldnames = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'LFP_prepro'};
for f = 1:length(fieldnames)
    Lfps.(fieldnames{f}) = Lfps.(fieldnames{f})(sesnum);
end