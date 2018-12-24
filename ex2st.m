function st = ex2st(ex, wnd, b)
% spike trains (trial x spike trains)
% INPUT: ex
%              wnd ... [start, end] of analysis window
%              b ... bin size (ms)
%
ex.Trials = ex.Trials(abs([ex.Trials.Reward]) > 0);
len_tr = length(ex.Trials);
ncol = round((1000/b)*(ex.fix.duration - wnd(2) - wnd(1)));
st = zeros(len_tr, ncol);
for i = 1:len_tr    
    t_strt = ex.Trials(i).Start - ex.Trials(i).TrialStart;
    
    begin = wnd(1);
    for c = 1:ncol
        st(i,c) = sum(ex.Trials(i).Spikes(ex.Trials(i).Spikes > t_strt(1) + begin ...
            & ex.Trials(i).Spikes <= t_strt(1) + begin + 0.001*b) - t_strt(1));
        begin = begin + 0.001*b;
    end
end