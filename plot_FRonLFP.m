function plot_FRonLFP
% check whether firing rate modulated by contrast stimulus affects stLFP
% ampltude
%
% Test whether the magnitude of the spike triggered LFP is depending on the
% spiking activity. We use the data recorded with a 2s stimulus, to avoid
% dominant slow fluctuations.
%
% 
% 04.04.18 Katsuhisa wrote it

close all

if nargin < 1; internal = 0; end

%% folder specifications
if mean(ismember('gpfs0', cd))==1
    main_dir = '/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/MagstLFP_vs_contrast/'; 
else
   main_dir = 'Z:\Katsuhisa\serotonin_project\LFP_project\MagstLFP_vs_contrast\';
end

%% load data given by 'check_FRonLFP.m'
load([main_dir 'Data/fr_lfp.mat'])

%% visualization
% stLFP in each session
lens = length([fr_lfp.session]);
wnd = fr_lfp.session(1).Results.wnd; % was 0.3
c = 0;
fignames = {'spikes vs contrast', 'LFP traces', 'STA', 'Power', 'Spectrogram'};
for i = 1:lens
    if isstruct(fr_lfp.session(i).Results)
        c = c + 1;
        
        % fr vs co in each unit
        figure(3);
        subplot(2,5,c)
%         nstm = length(fr_lfp.session(i).Results.stm.vals(2:end));
        me = fr_lfp.session(i).Results.spk_tu(2:end,1)';
        sem = fr_lfp.session(i).Results.spk_tu(2:end,2)'...
            ./fr_lfp.session(i).Results.sta.nspk(2:end);
        errorbar(fr_lfp.session(i).Results.stm.vals(2:end), me, sem, '-ok')
        
        if ismember(i, 5)
            ylabel('firing rate')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
        title([num2str(sum(fr_lfp.session(i).Results.ntr)) ' trials'])
        set(gca, 'XScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % LFP & stLFP traces
        col = [0 0 1; 0 1 0; 1 0 0];
        pl = nan(1,3);
        leg = cell(1, 3);
        for n = 1:3
            % LFP time course
            figure(4);
            subplot(2,5,c)
            lfpt = mean(fr_lfp.session(i).Results.lfp{end-3+n}, 1);
            lfpt = lfpt - mean(lfpt(fr_lfp.session(i).Results.ts <= 0));
            pl(n) = plot(fr_lfp.session(i).Results.ts, lfpt, 'color', col(n,:));
            leg{n} = [num2str(fr_lfp.session(i).Results.stm.vals(end-3+n))...
                ' (' num2str(fr_lfp.session(i).Results.sta.nspk(end-3+n)) ')' ];
            title([num2str(sum(fr_lfp.session(i).Results.ntr)) ' trials'])
            hold on;
            
            % STA 
            figure(5);
            subplot(2,5,c)
            plot(-wnd:0.001:wnd, fr_lfp.session(i).iResults.sta.mean(end-3+n, :), '-', 'color', col(n,:));
            hold on;
            title([num2str(sum(fr_lfp.session(i).iResults.ntr)) ' trials']) 
            
            % Power
            figure(6);
            subplot(2, 5 ,c)
            t_lfp = fr_lfp.session(i).Results.spectrogram.t{end-3+n};
            t_lfp = linspace(-0.2, 2, length(t_lfp)); 
            MeanS = mean(fr_lfp.session(i).Results.spectrogram.S{end-3+n}...
                (t_lfp > 1 & t_lfp <= 2, :), 1); 
            freq = fr_lfp.session(i).Results.spectrogram.f{end-3+n};
            plot(freq, 10*log10(MeanS), '-', 'color', col(n,:))
            hold on;
            
            % Spectrogram
            if n==3
                figure(7);
                subplot(2, 5, c)            
                spc = 10*log10(fr_lfp.session(i).Results.spectrogram.S{end-3+n});
                spc = spc - repmat(mean(spc(t_lfp < 0, :), 1), length(t_lfp), 1);
                imagesc(t_lfp, freq, spc')
                hold on;
                yy = get(gca, 'YLim');
                plot([0 0], yy, '--w')
                title([num2str(sum(fr_lfp.session(i).Results.ntr)) ' trials']) 
            end
        end
        for j = 1:2
            figure(j+3);
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            hold on;
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':k')
            ylim(yy)
            switch j
                case 1
                    xx = get(gca, 'XLim');
                    for l = 1:length(leg)
                        text(xx(1)+0.5*(xx(2)-xx(1)), yy(1)+0.1*(l-1)*(yy(2)-yy(1)), ...
                            leg{l}, 'color', col(l,:))
                    end
                    xlim([fr_lfp.session(i).Results.ts(1) fr_lfp.session(i).Results.ts(end)])
                case 2
                    xlim([-wnd wnd])
                    set(gca, 'XTick', [-wnd 0 wnd])
            end
        end
    end
end

% autosave
for i = 1:5
    h = figure(i+2);
    set(h, 'Name', fignames{i}, 'NumberTitle', 'off')
    savefig(h, [main_dir 'Figures/' fignames{i} num2str(internal) '.fig'])
end
disp('Figures saved!')
