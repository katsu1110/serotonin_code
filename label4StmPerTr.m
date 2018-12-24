function label_seq = label4StmPerTr(ex)
% assigns label of where the trial is positioned in the stimulus sequence:
% 0; fixation breaks
% 1, 2, 3...; the stimulus presentation number per trial
% written by Katsuhisa (09.08.18)

maxlab = ex.exp.StimPerTrial;
[~,idx] = sort(arrayfun(@(x) x.times.fpOn, ex.Trials));
oTrials = ex.Trials(idx);
len_tr = length(oTrials);
rewards = [oTrials.Reward];
label_seq = rewards;

% loop through all trials
try
    ct = 0;
    for i = 2:len_tr
        if rewards(i-1)==0 && rewards(i)==1
            ct = oTrials(i).Start(end) - oTrials(i).Start(1);
            label_seq(i) = 1;
        elseif rewards(i-1)==1 && rewards(i)==1
            ct = ct + oTrials(i).Start(end) - oTrials(i).Start(1);
            if label_seq(i-1) < maxlab && ct < 0.875*ex.fix.duration % consider blank
                label_seq(i) = label_seq(i-1) + 1;
            else
                ct = 0;            
                label_seq(i) = 1;
            end
        end
    end
catch
    return
end










