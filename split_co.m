function [ex0, ex2, idxs] = split_co(ex, es)
%%
% split ex-file with different contrast such that the effect size
% corresponds to 'es'
%

cos = unique([ex.Trials.co]);
lenc = length(cos);
spks = ones(1, lenc);
si = cell(1, lenc);
for c = 1:lenc
    [spkv, spkc] = getSpks(ex.Trials([ex.Trials.co]==cos(c)), [0.8 0]);
    spks(c) = sum(spkc);
    
    for i = 1:length(spkv)
        si{c}(i) = spktrain2psi(spkv{i}, [0.8 2], 50, 1000);
    end
end
oneidx = find(cos==1);
exes = spks/spks(oneidx);
[~, idx] = sort(abs(exes - es));
if idx(1)==oneidx
    idx = idx(2);
else
    idx = idx(1);
end
ex0 = ex;
ex0.Trials = ex.Trials([ex.Trials.co]==1);
ex0.si = si(cos==1);
ex2 = ex;
ex2.Trials = ex.Trials([ex.Trials.co]==cos(idx));
ex2.si = si(cos==cos(idx));
idxs = {find([ex.Trials.co]==1), find([ex.Trials.co]==cos(idx))};
