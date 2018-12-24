function [S,t,f] = lfpspecgram(exlfp, movingwin, params, k)
% plots the average stimulus spectogram
% in order to compare the stimulus induced changes, we 

[stimparam, vals] = getStimParam( exlfp );

% stimulus
trials = exlfp.Trials( [exlfp.Trials.(stimparam)] == vals(k));
[S, t, f] = computeSpectogram(trials , movingwin, params);

t = t+exlfp.Trials(1).LFP_prepro_time(1);

% if any(vals>1000) && vals(k) < 1000 
%     % blank
%     trials = exlfp.Trials( [exlfp.Trials.(stimparam)] == vals(vals>1000));
%     S_blank = computeSpectogram(trials , movingwin, params);
%     % normalize to a 'baseline' LFP
%     S = S-S_blank;
%     S(S<0)=0;
%     spec_helper([], [], exlfp, S, t, f, k);
% end

    spec_helper([], [], exlfp, S, t, f, k);


end


function [S,t, f] = computeSpectogram(trials, movingwin, params)
% compute the specogram given all the repeatadly recorded LFP

lfp = vertcat(trials.LFP_prepro);
lfp = lfp (mean(isnan(lfp), 2)==0,:);
params.err = 0;
[S,t,f] = mtspecgramc(lfp', movingwin, params );


end


function [] = spec_helper(~, ~, exlfp, S, t, f, k)
% plot the spectogram and name labels 

plot_matrix(S, t, f); % calls imagesc with logarithmically scaled power

colormap('jet')

[~, vals] = getStimParam( exlfp );
col = lines(length(vals));

title(['stim = ' num2str(vals(k))], 'Color', col(k,:));
xlabel('time [s]');
ylabel('frequency');

set(gca, 'FontSize', 8, 'XTick', [-0.2 0:0.2:t(end)], 'YTick', floor(f),...
    'CLim', [min(min(mag2db(S))) max(max(mag2db(S)))]);
xlim([t(1) t(end)]);

end