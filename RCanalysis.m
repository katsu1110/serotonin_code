% analysis script for reverse correlation data
clear
close all
clc

load('Z:\Katsuhisa\LFP_project\Data\rcdat.mat'); % dat struct including all fields from exinfo for rc data

figdir1 = fullfile('FiguresRC', 'SpikesandLFP');
mkdir(figdir1);

figdir2 = fullfile('FiguresRC', 'LFPandSTA');
mkdir(figdir2);
statime = 0.4;

for unit = 1:length(dat)
    
    try
       
        fprintf('\n\n==========================\n working on %1.1f \n', ...
            (dat(unit).id));
        %% ----------------------------------------------------- load ex files
        exbase = loadCluster(dat(unit).fname, 'ocul', dat(unit).ocul);
        exdrug = loadCluster(dat(unit).fname_drug, 'ocul', dat(unit).ocul);
        
        %% ---------------------------------------------------- lfp and spikes
%         h = figure('Name', dat(unit).figname);
%         plotSpikesAndLFP( exbase, 'xoff', 0.05, 'width', 0.4, 'fighandle', h);
%         plotSpikesAndLFP( exbase, 'xoff', 0.5, 'width', 0.4, 'fighandle', h);
%         
%         dat(unit).fig_spklfp = fullfile(figdir1, [dat(unit).figname '.fig']);
%         savefig(h, dat(unit).fig_spklfp, 'compact'); close(h);
        
        %% --------------------------- LFP in time and frequency domain and STA
        h = figure('Name', dat(unit).figname);
        
        % LFP in Time domain
        s1 = subplot(2,2,1);
        lfpTimeDomain(exbase); title('Baseline (k), 5HT (r)');
        s1.Children(1).Color = 'k'; s1.Children(1).ButtonDownFcn = [];
        s1.Legend.Location = 'northeast';
        s1.Legend.Visible ='off';
        
        
        h2 = figure; ax = axes; lfpTimeDomain(exdrug);
        ax.Children(1).Color = 'r';
        
        copyobj(ax.Children(1), s1); close(h2);
        
        % LFP in Frequency domain
        s2 = subplot(2,2,2);
        basep = lfpFreqDomain(exbase, [1 90]);
        s2.Children(1).Color = 'k'; s2.Children(1).ButtonDownFcn = [];
        s2.Children(2).FaceColor = 'k';
        
        h2 = figure; ax = axes; drugp=lfpFreqDomain(exdrug, [1 90]);
        ax.Children(1).Color = 'r';
        ax.Children(2).FaceColor = 'r';
        copyobj(ax.Children(1:2), s2); close(h2);

        f  = s2.Children(1).XData;
        y = drugp./basep;
        yyaxis right
        plot(f, y, ':'); 
        ylabel('Delta power (base-drug)')
        s2.YScale = 'log'; crossl;
        
        set(findobj(s1, 'Type', 'Text'), 'FontSize', 8);
        set(findobj(s2, 'Type', 'Text'), 'FontSize', 8);

        %%% save average power in dat
        dat(unit).alphapow_drug = mean(drugp(:,f>=7 & f<12)); % Alpha 7-12Hz
        dat(unit).alphapow = mean(basep(:,f>=7 & f<12));
        
        dat(unit).betapow_drug = mean(drugp(:,f>=12 & f<40)); %Beta 12-40Hz
        dat(unit).betapow = mean(basep(:,f>=12 & f<40));
        
        dat(unit).gammapow_drug = mean(drugp(:,f>=40 & f<90)); %Gamme 40-90Hz
        dat(unit).gammapow = mean(basep(:,f>=40 & f<90));
        
        
        title(['Avgerage pow (\alpha - \beta - \gamma)',...
            sprintf('\n  Baseline: %1.2f ', dat(unit).alphapow), ...
            sprintf('%1.2f ', dat(unit).betapow), ...
            sprintf('%1.2f \n', dat(unit).gammapow), ... 
            '  5HT : ' sprintf('%1.2f ', dat(unit).alphapow_drug), ...
            sprintf('%1.2f ', dat(unit).betapow_drug), ...
            sprintf('%1.2f  ', dat(unit).gammapow_drug)], ...
            'FontSize', 8, ...
            'Interpreter', 'tex');
        
        
        %%% spike triggered average LFP (STA)
        
        % Baseline
        s1 = subplot(2,2,3);
        [staLFP, ~, nspk]=spktriglfp(exbase, 'time', statime, 'plot'); hold on;
        [a1,t1]=getPeakAmplitude(staLFP, -statime:1/1000:statime);
        t1 = t1*1000;  % scale to ms
        
        s1.Children(1).ButtonDownFcn = [];
        set(findobj(gca, 'type', 'patch'), 'FaceColor', 'k');
        set(findobj(gca, 'type', 'line'), 'Color', 'k');
        
        % 5HT
        h2 = figure; ax = axes; %temporary objects
        [staLFP2, ~, nspk2]=spktriglfp(exdrug, 'time', statime, 'plot');
        [a2,t2]=getPeakAmplitude(staLFP2, -statime:1/1000:statime); 
        t2 = t2*1000; % scale to ms
        
        set(findobj(gca, 'type', 'patch'), 'FaceColor', 'r');
        set(findobj(gca, 'type', 'line'), 'Color', 'r');
        
        % superimpose the 5HT spike triggered LFP
        copyobj(ax.Children, s1);
        close(h2);
        
        set(findobj(gca, 'type', 'patch'), 'EdgeColor', 'w', ...
            'EdgeAlpha', 0, 'FaceAlpha', 0.4);
        crossl;
        title(sprintf('baseline(k) #spikes %1d \n5HT(r) #spikes %1d', ...
            nspk, nspk2));
        
         
        
        % show the STA in a smaller window
        s2= subplot(2,2,4);
        copyobj(s1.Children, s2); crossl;
        xlim([-0.05 0.05]); title('zoomed in STA');
        
        title(sprintf('baseline(k) amp %1.1fuV, t=%1.0fms \n5HT(r) amp %1.1fuV, t=%1.0fms', ...
            a1, t1, a2, t2));        
        
        
        set(findobj(s1, 'Type', 'Title'), 'FontSize', 8);
        set(findobj(s2, 'Type', 'Title'), 'FontSize', 8);
        
        dat(unit).fig_LFPandSTA = fullfile(figdir2, [dat(unit).figname '.fig']);
        savefig(h, dat(unit).fig_LFPandSTA, 'compact'); close(h);
        
        
        
        % assign infotmation to dat file
        dat(unit).stapeakamp = a1;
        dat(unit).stapeakt = t1;
        dat(unit).stapeakamp_drug = a2;
        dat(unit).stapeakt_drug = t2;
        
        dat(unit).stanspk = nspk;
        dat(unit).stanspk_drug = nspk2;
        
        
    catch
        fprintf('could not work on %s \n', dat(unit).figname);
        continue
    end
    
    %%
    close all;
    clearvars -except unit dat figdir1 figdir2 statime
end