function [fr_lfp] = check_FRonLFP
% check whether firing rate modulated by contrast stimulus affects stLFP
% ampltude
%
% Test whether the magnitude of the spike triggered LFP is depending on the
% spiking activity. We use the data recorded with a 2s stimulus, to avoid
% dominant slow fluctuations.
%
% 
% @CL 17.05.2017
% 04.04.18 Katsuhisa modified it

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

main_dir = [mypath '/Katsuhisa/serotonin_project/LFP_project/MagstLFP_vs_contrast/']; 
txt_fname = 'exlist.txt';   %<- file name containin the list of ex files 

%% initialize variables
% fig_pos = [788 552 959 420];
fr_lfp = struct('session', []);

%% load ex file names
fileID = fopen([main_dir txt_fname], 'r'); % read the text file
ex_fnames = textscan(fileID, '%s'); % scan the text and extract the file names
fclose(fileID);
ex_fnames = ex_fnames{1};
ex_fnames(5) = [];

%% compute the stLFP and its magnitude for each contrast for each unit 
% wnd = 0.3;
% t = -wnd:1/1000:wnd;
% idx = t>-0.015 & t<0.015;
out = cell(1, length(ex_fnames));
iout = cell(1, length(ex_fnames));
parfor i = 1:length(ex_fnames)
   % load file 
%    fname = strrep(ex_fnames{i}, 'c1', 'c0'); % cluster 0 <=> multiunit activity
   ex = loadCluster(ex_fnames{i}, 'loadlfp',1);
   
   % display some info
   disp('-------------------------------------------')
   disp(['fname: ' ex_fnames{i}])
   disp(['stm dur: ' num2str(ex.fix.stimDuration)])
   disp(['tf: ' num2str(ex.stim.vals.tf)])
   disp(['co: ' num2str(ex.exp.e1.range)])
   
   % preprocess
   period = {[0.5, 2]};
   
%    % filter LFP
%    ex = filterLFP(ex);
   ex.period = period;
      
   % LFP analysis
%    try
        out{i} = stmLFP(ex, 'LFP_prepro', {'all'}); 
        iout{i} = stmLFP(ex, 'iLFP_prepro', {'all'});
        disp(['Unit ' num2str(i) ' was processed.'])
%    catch
%        out{i} = nan;          
%        disp(['Unit ' num2str(i) ' had an error.'])
%    end
   
%    [stimdim, stimvals] = getStimParam(ex); %<- stimulus dimension and values
%    blank_i =  stimvals > 1;
%    ex_temp = ex;
%    
%    % stLFP in each stimulus type
%    lens = length(stimvals);
%    stalfp_m = nan(1, lens);
%    stalfp_t_m = nan(1, lens);
%    for s_i = 1:lens       
%        ex_temp.Trials = ex.Trials([ex.Trials.(stimdim)]==stimvals(s_i));
%        
%        av_stlfp = spktriglfp(ex_temp, 'plot', 'time', wnd);
%        % ??
%        av_stlfp = av_stlfp - mean(av_stlfp(t<-0.06));     
%        [stalfp_m(s_i), stalfp_t_m(s_i)] =  ...
%            getPeakAmplitude( av_stlfp(idx), t(idx) );       
%    end
%    
%    % store variables
%    fr_lfp.session(i).fname = ex_fnames{i};
%    fr_lfp.session(i).stmvals = stimvals(~blank_i);
%    fr_lfp.session(i).peak_stlfp = stalfp_m(~blank_i);
%    fr_lfp.session(i).t_peak_stlfp = stalfp_t_m(~blank_i);
   
%    % figures
%    h_fig = figure('Name', getFname(ex), 'Position', fig_pos);
%    
%    subplot(1,3,[1 2]);
%    
%    l = findobj(gca, 'Type', 'line');
%    col = lines(length(l));
%    for i = 1:length(l)
%        l(i).Color = col(i,:);
%    end   
%    legend('show');
%    crossl;
%    xlim([-0.3 0.3]);   
%    
%    s2 = subplot(1,3,3);
%    
%    % stlfp magnitude vs contrast
%    plot(stimvals(~blank_i), stalfp_m(~blank_i), 'ko-');  hold on;
%    if any(blank_i);   plot([eps 1], ones(1,2)*stalfp_m(blank_i), 'k--'); end
%    ylabel('stlfp magnitude (uV)')
%    
%    % time at peak/trough vs contrast
%    yyaxis(s2, 'right');
%    stalfp_t_m = stalfp_t_m*1000;
%    plot(stimvals(~blank_i), stalfp_t_m(~blank_i), 'bx-'); hold on;
%    if any(blank_i);   plot([eps 1], ones(1,2)*stalfp_t_m(blank_i), 'b--'); end
%    ylabel('stlfp peak t rel:stim onset (ms)');   
%    
%    k =  strfind(ex_fnames{i}, '\');
%    figname = ex_fnames{i}(k(end)+1:end);   
%     
%    axis square;
%    s2.XScale ='log';
%    s2.XLim = [min(stimvals(~blank_i)) 1];
%    s2.XLabel.String = 'contrast';
%    savefig(h_fig, fullfile([main_dir 'Figures'], ['stLFP_' figname '.fig']));
%    close(h_fig);
%       
% %    continue
%       
%    h_fig = figure('Name', ['tLFP_' getFname(ex)], 'Position', fig_pos);
%    lfpTimeDomain(ex);
%    savefig(h_fig, fullfile([main_dir 'Figures'], ['tLFP_' figname '.fig']));
%    close(h_fig);
% 
%    
%    h_fig = figure('Name', ['fLFP_' getFname(ex)], 'Position', fig_pos);
%    lfpFreqDomain(ex, [0 100]);
%    savefig(h_fig, fullfile([main_dir 'Figures'], ['fLFP_' figname '.fig']));
%    close(h_fig);

end
for i = 1:length(ex_fnames)
    fr_lfp.session(i).fname = ex_fnames{i};
    fr_lfp.session(i).Results = out{i};   
    fr_lfp.session(i).iResults = iout{i};
end

% autosave the output structure
save([main_dir 'Data/fr_lfp.mat'], 'fr_lfp', '-v7.3')
disp('Structure saved!')
