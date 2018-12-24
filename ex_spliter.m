function [exs, idx, val] = ex_spliter(ex, index, nsplit, stimparam)
% split ex file based on index
%

if nargin < 2; index = 'pupil'; end
if nargin < 3; nsplit = 2; end
if nargin < 4 
    [stimparam, stmvals] = getStimParam(ex);
else
    stmvals = unique([ex.Trials.(stimparam)]);
end

maxlab = max([ex.Trials.labelseq]);
maxidx = find([ex.Trials.labelseq]==maxlab, 1, 'last');
ntr = length(ex.Trials);
val = nan(1, ntr);
switch index
    case 'pupil'
        for i = 1:ntr
            if ex.Trials(i).labelseq==maxlab
                st = ex.Trials(i).Eye_prepro.stpos;
                en = ex.Trials(i).Eye_prepro.enpos;
                val(i - maxlab + 1:i) = ...
                    mean(ex.Trials(i).Eye_prepro.ps(st:en));
            end
        end    
        if maxidx < ntr
            st = ex.Trials(end).Eye_prepro.stpos;
            en = ex.Trials(end).Eye_prepro.enpos;
            val(maxidx+1:end) = mean(ex.Trials(end).Eye_prepro.ps(st:en));
        end
        idx = sortval(val, nsplit);
    case 'sc'
       % stimulus type
       lenv = length(stmvals);
       tr_stm = [ex.Trials.(stimparam)];
       idx = cell(1, nsplit);
       for i = 1:lenv
           % stm idx
           stmtr = find(tr_stm==stmvals(i));
           
           % spike counts
           [~, val] = getSpks(ex.Trials(stmtr), [0 0]);
           
           % sort already in each stimulus type
           sidx = sortval(val, nsplit);
           for n = 1:nsplit
                idx{n} = [idx{n}, stmtr(sidx{n})];
           end
       end   
end

exs = cell(1, nsplit);
for i = 1:nsplit
    exs{i} = ex;
    exs{i}.Trials = ex.Trials(idx{i});
end