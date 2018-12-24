%
% Test whether the magnitude of the spike triggered LFP is depending on the
% spiking activity. We use the data recorded with a 2s stimulus, to avoid
% dominant slow fluctuations.
%
%

% @CL 17.05.2017

close all
clc
%% folder specifications
prefix = 'Z:\Katsuhisa\serotonin_project\LFP_project\MagstLFP_vs_contrast\';
dat_dir = strcat(prefix, 'Data');     %<- directory for data
fig_dir = strcat(prefix, 'Figures');   %<- directory for figures
txt_fname = strcat(prefix, 'exlist.txt');   %<- file name containin the list of ex files 


%% initialize variables
fig_pos = [788 552 959 420];

%% load ex file names
fileID = fopen(txt_fname, 'r'); % read the text file
ex_fnames = textscan(fileID, '%s'); % scan the text and extract the file names
fclose(fileID);

ex_fnames = ex_fnames{1};


%% compute the stLFP and its magnitude for each contrast for each unit 
wnd =0.3;
t = -wnd:1/1000:wnd;

for unit_i = 1:length(ex_fnames)
    
   ex = loadCluster(ex_fnames{unit_i});
   
   
   [stimdim, stimvals] = getStimParam(ex); %<- stimulus dimension and values
   blank_i =  stimvals>1;
   ex_temp = ex;
   
   h_fig = figure('Name', getFname(ex), 'Position', fig_pos);
   
   subplot(1,3,[1 2]);
   for s_i = 1:length(stimvals)
       
       ex_temp.Trials = ex.Trials([ex.Trials.(stimdim)]==stimvals(s_i));
       
       av_stlfp = spktriglfp(ex_temp, 'plot', 'time', wnd);
       av_stlfp = av_stlfp- mean(av_stlfp(t<-0.06));
       
       idx = t>-0.015 & t<0.015;
       [stalfp_m(s_i), stalfp_t_m(s_i)] =  ...
           getPeakAmplitude( av_stlfp(idx), t(idx) );
       
   end
   
   
   l = findobj(gca, 'Type', 'line');
   col = winter(length(l));
   for i = 1:length(l)
       l(i).Color = col(i,:);
   end
   
   legend('show');
   crossl;
   xlim([-0.3 0.3]);
   
   
   s2 = subplot(1,3,3);
   
   % stlfp magnitude vs contrast
   plot(stimvals(~blank_i), stalfp_m(~blank_i), 'kx-');  hold on;
   if any(blank_i);   plot([eps 1], ones(1,2)*stalfp_m(blank_i), 'k--'); end
   ylabel('stlfp magnitude (uV)')

   
   % time at peak/trough vs contrast
   yyaxis(s2, 'right');
   stalfp_t_m = stalfp_t_m*1000;
   plot(stimvals(~blank_i), stalfp_t_m(~blank_i), 'bx-'); hold on;
   if any(blank_i);   plot([eps 1], ones(1,2)*stalfp_t_m(blank_i), 'b--'); end
   ylabel('stlfp peak t rel:stim onset (ms)');

   
   
   k =  strfind(ex_fnames{unit_i}, '\');
   figname = ex_fnames{unit_i}(k(end)+1:end);
   
    
   axis square;
   s2.XScale ='log';
   s2.XLim = [min(stimvals(~blank_i)) 1];
   s2.XLabel.String = 'contrast';
   savefig(h_fig, fullfile(fig_dir, ['stLFP_' figname '.fig']));
   close(h_fig);
   
   
   continue
      
   h_fig = figure('Name', ['tLFP_' getFname(ex)], 'Position', fig_pos);
   lfpTimeDomain(ex);
   savefig(h_fig, fullfile(fig_dir, ['tLFP_' figname '.fig']));
   close(h_fig);

   
   h_fig = figure('Name', ['fLFP_' getFname(ex)], 'Position', fig_pos);
   lfpFreqDomain(ex, [0 100]);
   savefig(h_fig, fullfile(fig_dir, ['fLFP_' figname '.fig']));
   close(h_fig);

end







