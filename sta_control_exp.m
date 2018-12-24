function stactrl = sta_control_exp(type)
%%
% visualize sta in the control experiment with different contrast in
% orientation reverse correlation stimuli
%

if nargin < 1; type = 'run'; end

% path ==========================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end

addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

folderpath = '/data/kaki/STAcontrol/';

switch type
    case 'run'
        % time relative to stimulus onset to filter and interpolate LFP
        t_off = -0.1;
        Fs = 1000;            % sampling frequency
        % stimulus presentation duration
        stimdur = 2;
        listings = dir([mypath folderpath]);
        listings(1:2) = [];
        lenf = length(listings);
        iout = cell(1, lenf);
        out = cell(1, lenf);
        cu = zeros(lenf, 2);
        wnd = 0.2;
        params = define_params(Fs, wnd, 1);
        movwin = [wnd, wnd/100];
        for i = 1:lenf
            % load data with LFP
            if contains(listings(i).name, 'lfp')
                % preprocess LFP
%                 ex = loadCluster([mypath folderpath listings(i).name], 'loadlfp', 1);
                ex = load([mypath folderpath listings(i).name]);
                ex = ex.ex;
                
                % load spikes
                exs = loadSpikes([mypath folderpath listings(i).name]);
                
                % combine
                oks = abs([ex.Trials.Reward]) > 0;
                ex.Trials = ex.Trials(oks);
                exs.Trials = exs.Trials(oks);
                lfpmat = [];
                stlfp = []; 
                ntr = length(ex.Trials);
                for n = 1:ntr
                    t_frame = ex.Trials(n).Start - ex.Trials(n).TrialStart; % time of frame onsets
                    t_lfp = ex.Trials(n).LFP_ts - t_frame(1) ; % time rel:stimulus onset
                    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
                    ex.Trials(n).LFP_prepro_time = time - t_frame(1);

                    % reduce the lfp signal to the period of stimulus presentation
                    ex.Trials(n).LFP_prepro = interp1(...
                        t_lfp, ex.Trials(n).LFP, ex.Trials(n).LFP_prepro_time);
                    ex.Trials(n).Spikes = exs.Trials(n).Spikes;
                    
                    % raw LFPs
                    lfpmat = [lfpmat; ex.Trials(n).LFP_prepro];
                    
                    % spectrogram
%                      [~, fs, ts, ps] = spectrogram_frange(ex.Trials(n).LFP_prepro, 90, Fs, [0 100]);
%                      if n==1
%                          S = ps;
%                      else
%                         S = S + ps;
%                      end
                     
                     % spike-triggered average LFP
                    spk = getSpks(ex.Trials(n), [0.2 0]);
                    try
                        temp = getSTA(ex.Trials(n).LFP_prepro, ex.Trials(n).LFP_prepro_time, ...
                                spk{1}, wnd, Fs);
                        stlfp = [stlfp; temp];
                    catch
                        disp('')
                    end
                end
                % spectrogram
                [S, ts, fs] = mtspecgramc(lfpmat', movwin, params);
                
                disp([listings(i).name ' loaded'])              

                % contrast, unit
                cu(i, :) = [ex.stim.vals.co_range, str2double(listings(i).name(6:7))];
%                 iout{i} = stmLFP(ex, 'iLFP_prepro', {'sta', 'spectrogram'});
%                 out{i} = stmLFP(ex, 'LFP_prepro', {'sta', 'spectrogram'});
                para.lfp = lfpmat;
                para.ts = ex.Trials(n).LFP_prepro_time;
                para.wnd = wnd;
                para.ntr = ntr;
                para.sta = stlfp;
                para.spectrogram = struct('f', fs, 't', ts, 'S', S');
                out{i} = para;
            end
        end
        out(cu(:, 1)==0) = [];
%         iout(cu(:, 1)==0) = [];
        cu(cu(:,1)==0, :) = [];
        stactrl.info = cu;
        stactrl.results = out;
%         stactrl.results0 = out;
%         stactrl.results1 = iout;
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/control/stactrl.mat'])
        disp('saved!')
    case 'plot'
        % load struct
        load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/control/stactrl.mat'], 'stactrl')
        
        % unit
        uniuni = unique(stactrl.info(:,2));
        lenu = length(uniuni);
        
        % visualize
        close all;
        figure;
        cols = [1 0 0; 0 1 0; 0 0 1; 0 1 1];
        for i = 1:length(uniuni)
            corange = sort(stactrl.info(stactrl.info(:,2)==uniuni(i), 1), 'descend');
            nspk = zeros(1, length(corange));
            for c = 1:length(corange)
                idx = find(stactrl.info(:,1)==corange(c) & stactrl.info(:,2)==uniuni(i));
                idx = idx(1);
                
                % LFP traces
                subplot(5, lenu, i)
                me = nanmean(stactrl.results{idx}.lfp, 1);
                ts = stactrl.results{idx}.ts;
                plot(ts, me, '-', 'color', cols(c, :))
                hold on;
                
                % Spectrogram
                if c==1
                    subplot(5, lenu, i + lenu)
                    tp = stactrl.results{idx}.spectrogram.t;
                    tp = linspace(-0.1, 2, length(tp));
                    f = stactrl.results{idx}.spectrogram.f;
                    S = 10*log10(stactrl.results{idx}.spectrogram.S);
                    S = S - repmat(mean(S(:, tp < 0), 2), 1, length(tp));
                    imagesc(tp, f, S)
                    hold on;
                end
                                
                % power
                subplot(5, lenu, i + 2*lenu)
                S = 10*log10(stactrl.results{idx}.spectrogram.S);
%                 S = S - repmat(mean(S(:, tp < 0), 2), 1, length(tp));
                plot(f, mean(S(:, tp > 0.5), 2), '-', 'color', cols(c, :))
                hold on;
                
                subplot(5, lenu, i + 3*lenu)
%                 S = 10*log10(stactrl.results{idx}.spectrogram.S);
                S = S - repmat(mean(S(:, tp < 0), 2), 1, length(tp));
                plot(f, mean(S(:, tp > 0.5), 2), '-', 'color', cols(c, :))
                hold on;
                
                % STA
                subplot(5, lenu, i + 4*lenu)
                me = nanmean(stactrl.results{idx}.sta, 1);
%                 sd = nanstd(stactrl.results{idx}.sta, [], 1);
                nspk(c) = size(stactrl.results{idx}.sta, 1);
%                 sem = sd/sqrt(nspk(c));
                t = linspace(-stactrl.results{idx}.wnd, stactrl.results{idx}.wnd, length(me));
%                 fill_between(t, me - sem, me + sem, cols(c, :), 0.5);
                hold on;
                plot(t, me, '-', 'color', cols(c, :))
                hold on;          
            end
            
            % format
            subplot(5, lenu, i)
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':k')        
            axis([ts(1) ts(end) yy])
            title(['unit ' num2str(uniuni(i))])
            if i==1
                xlabel('time from stimulus onset (sec)')
                ylabel('LFP')
            end
            set(gca, 'box', 'off', 'tickdir', 'out')
            
            subplot(5, lenu, i + lenu)
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':w')        
            axis([tp(1) tp(end) yy])
            if i==1
                xlabel('time from stimulus onset (sec)')
                ylabel('frequency (Hz)')
            end
            set(gca, 'box', 'off', 'tickdir', 'out')
            
            subplot(5, lenu, i + 2*lenu)
            if i==1
                ylabel('power')
                xlabel('frequency (Hz)')
            end
%             set(gca, 'YScale', 'log')
            set(gca, 'box', 'off', 'tickdir', 'out')
            
            subplot(5, lenu, i + 3*lenu)
            if i==1
                ylabel('\Delta power')
                xlabel('frequency (Hz)')
            end
%             set(gca, 'YScale', 'log')
            set(gca, 'box', 'off', 'tickdir', 'out')
            
            subplot(5, lenu, i + 4*lenu)
            xx = get(gca, 'XLim');
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':k')
            for c = 1:length(corange)
                text(xx(1)+0.6*(xx(2)-xx(1)), yy(1)+(0.95 - 0.1*(c-1))*(yy(2)-yy(1)), ...
                    [num2str(corange(c)) ' (' num2str(nspk(c)) ')'], 'color', cols(c, :))
            end
            axis([xx yy])
            xlim([-0.03 0.03])
            if i==1
                xlabel('time from spike (sec)')
                ylabel('LFP')
            end
            set(gca, 'box', 'off', 'tickdir', 'out')
        end        
        
%     case 'plot'
%         % load struct
%         load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/control/stactrl.mat'], 'stactrl')
%         
%         % unit
%         uniuni = unique(stactrl.info(:,2));
%         lenu = length(uniuni);
%         
%         % visualize
%         close all;
%         figure;
%         cols = [1 0 0; 0 1 0; 0 0 1; 0 1 1];
%         for i = 1:length(uniuni)
%             corange = sort(stactrl.info(stactrl.info(:,2)==uniuni(i), 1), 'descend');
%             nspk = zeros(1, length(corange));
%             for c = 1:length(corange)
%                 idx = find(stactrl.info(:,1)==corange(c) & stactrl.info(:,2)==uniuni(i));
%                 idx = idx(1);
%                 
%                 % LFP traces
%                 subplot(4, lenu, i)
%                 me = nanmean(stactrl.results0{idx}.lfp{end}, 1);
%                 ts = stactrl.results0{idx}.ts;
%                 plot(ts, me, '-', 'color', cols(c, :))
%                 hold on;
%                 
%                 % Spectrogram
%                 if c==1
%                     subplot(4, lenu, i + lenu)
%                     tp = stactrl.results0{idx}.spectrogram.t{end};
%                     tp = linspace(-0.2, 2, length(tp));
%                     f = stactrl.results0{idx}.spectrogram.f{end};
%                     S = 10*log10(stactrl.results0{idx}.spectrogram.S{end})';
%                     S = S - repmat(mean(S(:, tp < 0), 2), 1, length(tp));
%                     imagesc(tp, f, S)
%                     hold on;
%                 end
%                                 
%                 % power
%                 subplot(4, lenu, i + 2*lenu)
%                 S = 10*log10(stactrl.results0{idx}.spectrogram.S{end})';
% %                 S = S - repmat(mean(S(:, tp < 0), 2), 1, length(tp));
%                 plot(f, mean(S(:, tp > 0.5), 2), '-', 'color', cols(c, :))
%                 hold on;
%                 
%                 % STA
%                 subplot(4, lenu, i + 3*lenu)
%                 me = stactrl.results1{idx}.sta.mean;
%                 nspk(c) = stactrl.results1{idx}.sta.nspk;
%                 sem = stactrl.results1{idx}.sta.sd/sqrt(nspk(c));
%                 t = linspace(-0.05, 0.05, length(me));
%                 fill_between(t, me - sem, me + sem, cols(c, :), 0.5);
%                 hold on;
%                 plot(t, me, '-', 'color', cols(c, :))
%                 hold on;          
%             end
%             
%             % format
%             subplot(4, lenu, i)
%             yy = get(gca, 'YLim');
%             plot([0 0], yy, ':k')        
%             axis([ts(1) ts(end) yy])
%             title(['unit ' num2str(uniuni(i))])
%             if i==1
%                 xlabel('time from stimulus onset (sec)')
%                 ylabel('LFP')
%             end
%             set(gca, 'box', 'off', 'tickdir', 'out')
%             
%             subplot(4, lenu, i + lenu)
%             yy = get(gca, 'YLim');
%             plot([0 0], yy, ':w')        
%             axis([tp(1) tp(end) yy])
%             if i==1
%                 xlabel('time from stimulus onset (sec)')
%                 ylabel('frequency (Hz)')
%             end
%             set(gca, 'box', 'off', 'tickdir', 'out')
%             
%             subplot(4, lenu, i + 2*lenu)
%             if i==1
%                 ylabel('power')
%                 xlabel('frequency (Hz)')
%             end
% %             set(gca, 'YScale', 'log')
%             set(gca, 'box', 'off', 'tickdir', 'out')
%             
%             subplot(4, lenu, i + 3*lenu)
%             xx = get(gca, 'XLim');
%             yy = get(gca, 'YLim');
%             plot([0 0], yy, ':k')
%             for c = 1:length(corange)
%                 text(xx(1)+0.6*(xx(2)-xx(1)), yy(1)+(0.95 - 0.1*(c-1))*(yy(2)-yy(1)), ...
%                     [num2str(corange(c)) ' (' num2str(nspk(c)) ')'], 'color', cols(c, :))
%             end
%             axis([xx yy])
%             if i==1
%                 xlabel('time from spike (sec)')
%                 ylabel('LFP')
%             end
%             set(gca, 'box', 'off', 'tickdir', 'out')
%         end        
end

% subfunction
function ex = loadSpikes(fname)
fname = strrep(fname, 'lfp', 'c1');
ex = load(fname);
ex = ex.ex;

