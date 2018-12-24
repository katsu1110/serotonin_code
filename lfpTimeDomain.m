function lfp_avg = lfpTimeDomain(exlfp)
% stimulus triggered average

[ stimparam, vals] = getStimParam( exlfp );
col = lines(length(vals));

col = flipud(col);

ts = exlfp.time;


% exlfp.Trials = exlfp.Trials([exlfp.Trials.phase]==4);

%%%
for i = 1:length(vals)
    
    trials = exlfp.Trials([exlfp.Trials.(stimparam)] == vals(i));

    %%% averaged trials
    lfp_avg(i,:) = nanmean(vertcat(trials.LFP_prepro), 1);
        
    plot(ts, lfp_avg(i,:), 'Color', col(i, :),...
        'DisplayName', num2str(vals(i)), ...
        'ButtonDownFcn',  ...
        {@PlotSingleTrials, trials, vals(i), ts}); hold on
end

% legend('show', 'Location', 'southoutside', 'Orientation', 'horizontal');
xlim([ts(1) ts(end)]);
xlabel('time (s)'); ylabel('avg LFP (\muV)');
crossl; 
end




function PlotSingleTrials(source, ~,trials, stim, ts)
% Callback function for single averaged signals to show the individual
% trials

figure('Name', ['All Time Domain Trials Stim=' num2str(stim)]);
off = 0;

for i= 1:length(trials)
    
    
    plot([ts(1) ts(end)], [off, off], 'k:'); hold on;
    
    % filtered signal
    if trials(i).FixCounter == 1
        lw = 1;
    else
        lw = 0.5;
    end
    
    p = plot(trials(i).LFP_prepro_time, trials(i).LFP_prepro+off, ...
        'Color', source.Color, 'LineWidth', lw); hold on;
        
    p.Color = [p.Color 0.1]; % make the lines transparent
%     off = off + max(trials(i).LFP_prepro)*2;
end
    
plot([0 0], get(gca, 'YLim'), 'k--');
title(sprintf('lfp in trials of stimulus %1.2f \n %d repeats', stim, length(trials)));
xlabel('time');
ylabel('stacked LFP');
xlim([ts(1) ts(end)]);
end



