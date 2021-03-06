function [ex0, ex2] = zscore_pairLFP(ex0, ex2, lfpfield)

% zscoring
lfpall = []; 
for i = 1:length(ex0.Trials)
    lfpall = [lfpall, ex0.Trials(i).(lfpfield)];
end
for i = 1:length(ex2.Trials)
    lfpall = [lfpall, ex2.Trials(i).(lfpfield)];
end
me = nanmean(lfpall);
sd = nanstd(lfpall);
for i = 1:length(ex0.Trials)
    ex0.Trials(i).(lfpfield) = (ex0.Trials(i).(lfpfield) - me)/sd;
end
for i = 1:length(ex2.Trials)
    ex2.Trials(i).(lfpfield) = (ex2.Trials(i).(lfpfield) - me)/sd;
end
