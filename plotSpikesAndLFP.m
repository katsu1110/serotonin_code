function figh = plotSpikesAndLFP( ex, varargin )
% plot the spike train together with two bandpass filtered versions of the 
% LFP for a subset of rewarded trials
%
% additional arguments are
% ntrials - the first x trials are shown in the graph
% alpharg - the first frequency band shown in the plot
% gammarg - the second frequency band shown in the plot 
%
% @CL 10.05.2017


%% define variables

Fs = 1000; % Sampling frequency
% frequency range for the  two signals shown in the graph
alpharg = [1 4]; % alpha range
gammarg = [20 150]; % gamma range

% trials shown
ntrials = 10;


xoff = 0.05; % axis offset on the left side
w = 1-2*xoff; % axis width

%% parse additional input
j = 1;

while j<= length(varargin)
    switch varargin{j}
        case 'ntrials'
            ntrials = varargin{j+1};
        case 'alpharg'
            alpharg = varargin{j+1};
        case 'gammarg'
            gammarg = varargin{j+1};
        case 'xoff' 
            xoff = varargin{j+1};
        case 'width'
            w = varargin{j+1};
        case 'fighandle'
            figh = varargin{j+1};
    end
    j = j+2;
end
     
%% design filter

[b_alpha, a_alpha] = butter(2, alpharg/(Fs/2), 'bandpass');
[b_gamma, a_gamma] = butter(2, gammarg/(Fs/2), 'bandpass');

%% start building the figure
if ~exist('figh', 'var')
    figh = figure('Position', [680 60 560 924], 'Name', getFname(ex));
else
    figh.Position = [680 60 560 924];
end

% axes are equally spaced for each trial
yoff = 0.05;
h = (1-2*yoff)/(ntrials*3 + (ntrials-1) ); % axes height
for trial = 1:ntrials
    
    % extract spikes, lfp and trial time stamps
    spks = num2cell(ex.Trials(trial).Spikes);
    lfp = ex.Trials(trial).LFP_prepro;
    t = ex.Trials(trial).LFP_prepro_time;
    
    % plot spike train
    s(trial*3-2) = axes('Position', [xoff yoff+h w h/2]); 
    hold on;
    cellfun(@(x) plot([x, x], [0.2 0.7], 'r', 'LineWidth', 0.5), spks, 'UniformOutput', 0); 
    axis off; ylim([0 1]);
    
    % plot 'alpha' oscillation
    s(trial*3-1) = axes('Position', [xoff yoff+h w h*1.5]);  
    plot(t, filtfilt(b_alpha, a_alpha, lfp), 'k'); 
    axis off; crossl;

    % plot 'gamma' oscillation
    s(trial*3) = axes('Position', [xoff yoff w h*1.5]); 
    plot(t, filtfilt(b_gamma, a_gamma, lfp), 'b'); 
    axis off; crossl;
    
    yoff = yoff+h*4;        
end

set(s, 'XLim', [-0.15, t(end)], 'Clipping', 'on');

linkaxes(s(2:3:end), 'y'); % link axes with corresponding LFP frequency bands
linkaxes(s(3:3:end), 'y');

% adjust the first two axes to show the scaling
set(s(2:3), 'XTick', [],'Box', 'off', ...
    'Color','none', 'XColor','none','YColor','k', ...
     'TickDir', 'out');
set(s(3), 'YAxisLocation','right', 'XColor','k'); 

s(2).YTick = s(2).YLim; 
s(3).YTick = s(3).YLim; s(3).XTick = sort([0, s(3).XLim]);

axes(s(3)); axis on; xlabel('time (s)');
axes(s(2)); axis on;

%% information regarding the shown graphics
axes('Position', [xoff yoff+h/2 w h/2]); axis off;


text(0,0, sprintf('RED spikes   BLACK LFP %1d-%1d Hz band  \nBLUE LFP %1d-%1d Hz band in uV', ...
    alpharg, gammarg));

set(findobj(gcf, 'Type', 'Text'), 'FontSize', 8);

end

