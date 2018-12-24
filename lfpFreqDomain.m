function pow_avg = lfpFreqDomain(exlfp, frange)
% stimulus triggered average within a particular frequency range frange

[ stimparam, vals] = getStimParam( exlfp );
col = lines(length(vals));
col = flipud(col);


% exlfp.Trials = exlfp.Trials([exlfp.Trials.phase]==4);


%%% frequency domain offset
fidx = exlfp.Trials(1).FREQ>= frange(1) & exlfp.Trials(1).FREQ <= frange(2);
f = exlfp.Trials(1).FREQ(fidx);

%%%
for i = 1:length(vals)
    
    trials = exlfp.Trials([exlfp.Trials.(stimparam)] == vals(i));
    
    % average power across the same stimulus condition 
    pow_avg(i,:) = nanmean(horzcat(trials.POW),2);
    pow_sem(i,:) = nanstd(horzcat(trials.POW),0, 2)/sqrt(length(trials));
    
    
    fill([f; flipud(f)], [pow_avg(i,fidx)+pow_sem(i,fidx),...
        fliplr(pow_avg(i,fidx)-pow_sem(i,fidx))], col(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    
    plot(f, pow_avg(i,fidx), 'Color', col(i, :),...
        'DisplayName', num2str(vals(i)),...
        'ButtonDownFcn', ...
        {@PlotSingleTrials, f, fidx, trials, vals(i), col(i,:)}); hold on
    
end
xlim([f(1) f(end)]);
xlabel('Frequency (Hz)'); ylabel('avg Power +/- SEM (\muV^2)');
set(gca, 'YScale', 'log');
pow_avg = pow_avg(:,fidx); 
end


function PlotSingleTrials(source, ~,f, fidx, trials, stim, col)
% Callback function for single averaged signals to show the individual
% trials

figure('Name', ['All Frequency Domain Trials Stim=' num2str(stim)]);
off = 0;
for i= 1:length(trials)
    
    plot(f, log(trials(i).POW(fidx))+off, 'Color', source.Color); hold on;
    
%     off = off + max(log(trials(i).POW(fidx)))*2;
end

title(['trial for stimulus = ' num2str(stim)])
xlabel('Frequency');
ylabel('stacked LFP');

end



