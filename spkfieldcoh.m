function [C, f] = spkfieldcoh(ex, params)

[ stimparam, vals] = getStimParam( ex );

col = lines(length(vals));

for i = 1:length(vals)
    
    idx = [ex.Trials.(stimparam)] == vals(i);
    lfp = vertcat(ex.Trials(idx).LFP_prepro);
    spk = getSpks(ex.Trials(idx));
    
    [C(i,:),phi,S12,S1,S2,f(i,:)]= coherencycpt(lfp',spk,params);
    plot(f(i,:), C(i,:), 'Color', col(i,:), 'DisplayName', num2str(vals(i))); hold on;
end

ylabel('Coherence');
xlabel('Frequency');
legend('show', 'Location', 'EastOutside')

end


function spk_out = getSpks(trial)

for i = 1:length(trial)
% spikes within the stimulus presentation time
t_strt = trial(i).Start - trial(i).TrialStart;
t_end = t_strt(end)+mean(diff(t_strt));

spk = trial(i).Spikes( trial(i).Spikes >= t_strt(1) & ...
    trial(i).Spikes <= t_end) - t_strt(1); 
spk = round(spk*1000)/1000;

spk_out(i).spk = spk;
end
end

