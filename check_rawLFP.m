function check_rawLFP
%%
% check raw LFP data
%

% list of data ======================
datalist = cell(1, 4);
% Kaki, recent, rig2
datalist{1} = 'Z:\data\kaki\STAcontrol\ka_0014_c1_sortKK_12.11.grating.ORxRC.mat';
% Kaki, old, rig1
datalist{2} = 'Z:\data\kaki\0131\ka_0131_c1_sortLH_12.30.grating.ORxRC.mat';
% Kiwi, recent, rig2
datalist{3} = 'Z:\data\kiwi\sorted\0329\ki_0329_c1_sortLH_16.59.rds.DXxDCxDC2xDX2.mat';
% Mango, old, rig1
datalist{4} = 'Z:\data\mango\0106\ma_0106_c1_sortHN_3.06PM.grating.ORxRC.mat';

titlelab = {'Rig2 - recent (kaki)', 'Rig1 -old (Kaki)', 'Rig2 -recent (kiwi)', 'Rig1 -old (mango)'};

% LFP traces ======================
close all;
% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.1;
Fs = 1000;            % sampling frequency
% stimulus presentation duration
stimdur = 2;
for d = 1:length(datalist)
    % load LFP data
    % preprocess LFP
    ex = loadCluster(datalist{d}, 'loadlfp', 1);
%     fname = strrep(datalist{d}, '_c1_', '_lfp_');
%     ex = load(fname);
%     ex = ex.ex;
    
    % load spike data
    ex2 = load(datalist{d});
    
    % plot LFP traces
    Trials = ex.Trials(abs([ex.Trials.Reward]) > 0);
    sTrials = ex2.ex.Trials(abs([ex.Trials.Reward]) > 0);
    ntr = length(Trials);
    lfpmat = [];
    stlfp = []; 
    wnd = 0.2;
    for n = 1:ntr        
%         % LFP data =======================
%         t_frame = Trials(n).Start - Trials(n).TrialStart; % time of frame onsets
%         t_lfp = Trials(n).LFP_ts - t_frame(1) ; % time rel:stimulus onset
%         time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
%         Trials(n).LFP_prepro_time = time - t_frame(1);
% 
%         % reduce the lfp signal to the period of stimulus presentation
%         Trials(n).LFP_prepro = interp1(...
%             t_lfp, Trials(n).LFP, Trials(n).LFP_prepro_time);
        lfpmat = [lfpmat; Trials(n).LFP_prepro];
        
        % spectrogram
         [~, fs, ts, ps] = spectrogram_frange(Trials(n).LFP_prepro, 90, Fs, [0 100]);
         if n==1
             S = ps;
         else
            S = S + ps;
         end
          
        % spike-triggered average LFP
%         spk = getSpks(sTrials(n), [0.5 0]);
        spk = getSpks(Trials(n), [0.5 0]);
        stlfp = [stlfp; getSTA(Trials(n).LFP_prepro, Trials(n).LFP_prepro_time, ...
                spk{1}, wnd, Fs)];
        
        % plot
        figure(1);
        subplot(2, 2, d)
        p = plot(Trials(n).LFP_prepro_time, Trials(n).LFP_prepro, '-', ...
            'color', 0.1*[1 1 1], 'linewidth', 0.25);
        p.Color(4) = 0.1;
        hold on;
    end

    % mean
    tl = Trials(end).LFP_prepro_time;
    plot(tl, mean(lfpmat, 1), '-', ...
            'color', 'r', 'linewidth', 3)
    xlabel('time after stimulus onset (sec)')
    ylabel('LFP (uV)')
    title(titlelab{d})
    yy = get(gca, 'YLim');
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    xlim([Trials(end).LFP_prepro_time(1), ...
        Trials(end).LFP_prepro_time(end)])
    ylim(yy)
    set(gca, 'box', 'off', 'tickdir', 'out')
    
    % power spectra
    S = S/ntr;
    figure(2);    
    subplot(2, 2, d)
    ts = linspace(t_off, stimdur, length(ts));
    plot(fs, mean(S(:, 0.5 < ts & ts < 2), 2), '-k')
    hold on;
    set(gca, 'YScale', 'log')
    xlabel('frequency (Hz)')
    ylabel('power')
    title(titlelab{d})    
    set(gca, 'box', 'off', 'tickdir', 'out')
    
    % Spectrogram
    figure(3);
    subplot(2, 2, d)
    imagesc(ts, fs, 10*log10(S))
    hold on;
    yy = get(gca, 'YLim');
    plot([0 0], yy, '--w')
    xlabel('time (sec)')
    ylabel('frequency (Hz)')
    title(titlelab{d})    
    set(gca, 'box', 'off', 'tickdir', 'out')
    
    % delta Spectrogram
    figure(4);
    subplot(2, 2, d)
    Sres = 10*(log10(S) - repmat(mean(log10(S(:, ts < 0)), 2), 1, length(ts)));
    imagesc(ts, fs, Sres)
    hold on;
    yy = get(gca, 'YLim');
    plot([0 0], yy, '--w')
    xlabel('time (sec)')
    ylabel('frequency (Hz)')
    title(titlelab{d})    
    set(gca, 'box', 'off', 'tickdir', 'out')
    
    % STA
    figure(5);
    subplot(2, 2, d)
    t_sta = linspace(-wnd, wnd, size(stlfp, 2));
    me_sta = mean(stlfp, 1);
    nspk = size(stlfp, 1);
    sem_sta = std(stlfp, [], 1)/sqrt(nspk);
    fill_between(t_sta, me_sta - sem_sta, me_sta + sem_sta, ...
        [0 0 0], 0.5);
    hold on;
    plot(t_sta, me_sta, '-k', 'linewidth', 2)
    yy = get(gca, 'YLim');
    plot([0 0], yy, '--k')
    xlim([-wnd wnd])
    ylim(yy)
    text(0.05, yy(1)+0.95*(yy(2)-yy(1)), [num2str(nspk) ' spikes'])
    xlabel('time (sec)')
    ylabel('LFP (uV)')
    title(titlelab{d})    
    set(gca, 'box', 'off', 'tickdir', 'out')
    
end
figure(1);
set(gcf, 'Name', 'raw LFP traces (uV)', 'NumberTitle', 'off')
figure(2);
set(gcf, 'Name', 'spectra power', 'NumberTitle', 'off')
figure(3);
set(gcf, 'Name', 'Spectrogram', 'NumberTitle', 'off')
figure(4);
set(gcf, 'Name', 'Spectrogram (baseline corrected)', 'NumberTitle', 'off')
figure(5);
set(gcf, 'Name', 'spike triggered average LFP', 'NumberTitle', 'off')