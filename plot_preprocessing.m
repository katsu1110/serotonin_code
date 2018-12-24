function plot_preprocessing(ex)
% visualize the effects of preprocessing
% - stimulus triggered LFP
% - power of it
% - spike triggered LFP
% - power of it

close all;
h = figure;
ncol = length(ex.Trials(1).preproLFP);
N = length(ex.Trials);
disp(['The number of trials: ' num2str(N)])

% visualize
for c = 1:ncol    
    for n = 1:N
        % stimulus triggered LFP
        subplot(2, ncol, c)
        plot(ex.time, ex.Trials(n).preproLFP(c).lfp)
        hold on;
                
        % power
        subplot(2, ncol, c + ncol)
        plot(ex.Trials(n).preproLFP(c).FREQ, ex.Trials(n).preproLFP(c).POW)        
        hold on;
    end
    subplot(2, ncol, c)
    switch c
        case 1
            title('no preprocessing')
        case 2
            title('+ baseline correction')
        case 3
            title('+bandpass filtering')
        case 4
            title('+ detrending')
        case 5
            title('+ line noise removal')
        case 6
            title('+ notch filter')
    end    
    [~, zerot] = min(abs(ex.time));
    yy = get(gca, 'YLim');
    plot(zerot*[1 1], yy, '--k')
    xlim([ex.time(1) ex.time(end)])
    ylim(yy)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2, ncol, c + ncol)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    set(gca, 'YScale', 'log')
end