function ex = frequAnalysis( ex, varargin )
% filters and interpolates lfp traces and returns their the ex file also
% containing the filtered, interpolated and power spectra
%
% optional input variables are:
% fs : sampling frequency
% notchf : notch filtered frequency
% notchord : notch filter order
% bpf:   band pass filtered frequency
% bpord: band pass filter order
% nw    : time bandwidth product
%
% @CL 
%


%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% stimulus per trials
if ex.fix.duration == 2
    label_seq = label4stmPerTr(ex, 1);
else
    label_seq = label4stmPerTr(ex, 4);
end

% notch filter variables
Fs = 1000;              % sampling frequency
% notchf = [49 51];     % notch filter frequency1
% notchf2 = [99 101];
% notchf3 = [149 151];
% notchf4 = [33 35];
% notchf5 = [66 68];
% notchord = 2;           % filter order

% bpf = [0.5 100];     % bandpass filter cutoff frequency
bpord = 2;            % filter order

% % power spectrum in Chalk et al. time halfbandwidth was 3 and k was 5
% nw = 2;         % time half bandwidth product
% 
% nfft = 1024; % number of frequency samples - I think

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.2;

% % parameters for chronux toolbox
% params = define_params;

% TODO: how to handle files with missing LFP entries (such as ma0014,
% 5.44pm)
if any(cellfun(@isempty, {ex.Trials.LFP}))
%     warning('TODO: in frequAnalysis line 44. How to proceed with empty LFP entries.');
    label_seq = label_seq(~cellfun(@isempty, {ex.Trials.LFP}));
    ex.Trials = ex.Trials(~cellfun(@isempty, {ex.Trials.LFP}));
end

% %% parse input
% k = 1;
% while k<=length(varargin)
%     switch varargin{k}
%         case 'fs'
%             Fs = varargin{k+1};
%         case 'notchf'
%             notchf = varargin{k+1};
%         case 'notchord'
%             notchord = varargin{k+1};
%         case 'bpf'
%             bpf = varargin{k+1};
%         case 'bpord'
%             bpord = varargin{k+1};
%         case 'nw'
%             nw = varargin{k+1};
%     end
%     k=k+2;
% end
% clearvars k;

%% generate filters
% define notch filter
% [b_notch,a_notch] = butter(notchord, notchf/(Fs/2), 'stop' );
% [b_notch2,a_notch2] = butter(bpord, notchf2/(Fs/2), 'stop');
% [b_notch3,a_notch3] = butter(bpord, notchf3/(Fs/2), 'stop');
d1 = designfilt('bandstopiir','FilterOrder', bpord, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',Fs);
ltf = floor(ex.stim.vals.tf);
utf = ceil(ex.stim.vals.tf);
d2 = designfilt('bandstopiir','FilterOrder', bpord, ...
               'HalfPowerFrequency1', ltf, 'HalfPowerFrequency2', utf, ...
               'DesignMethod','butter','SampleRate',Fs);
d3 = designfilt('bandstopiir','FilterOrder', bpord, ...
               'HalfPowerFrequency1',99,'HalfPowerFrequency2',101, ...
               'DesignMethod','butter','SampleRate',Fs);
d4 = designfilt('bandstopiir','FilterOrder', bpord, ...
               'HalfPowerFrequency1',149,'HalfPowerFrequency2',151, ...
               'DesignMethod','butter','SampleRate',Fs);
           
% % define bandpass filter
% [b_bpass,a_bpass] = butter(bpord, bpf/(Fs/2), 'bandpass');
% [b_bpass_notch,a_bpass_notch] = butter(bpord/2, notchf/(Fs/2), 'bandpass');

%% perform functions on each trial lfp
ex.Trials = ex.Trials(label_seq > 0);
N = length(ex.Trials);
for n = 1:N    
    % timing information
    t_frame = ex.Trials(n).Start - ex.Trials(n).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(n).LFP_ts - t_frame(1) ; % time rel:stimulus onset
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    time = time - t_frame(1);

    % bandpass filter
%     filt1 = filtfilt(b_bpass, a_bpass, ex.Trials(ind).LFP);        
        filt1 = ex.Trials(n).LFP;
        
    %% notch filter - filters line noise    
    % apply notch filter
%     filt2 = filtfilt(b_notch, a_notch, filt1);
%     filt2 = filtfilt(b_notch2, a_notch2, filt2);
%     filt2 = filtfilt(b_notch3, a_notch3, filt2);
%     filt2 = filtfilt(d, filt1);
    filt2 = filtfilt(d1, filt1);
    filt2 = filtfilt(d2, filt2);
    filt2 = filtfilt(d3, filt2);
    filt2 = filtfilt(d4, filt2);
    
%     % remove 50Hz line noise with regression
%     filt2 = rmlinesc(filt2,params,0.05/N,0,50);

%     filt2 = ex.Trials(ind).LFP; % for debugging only
    % detrending
%     filt2 = locdetrend(filt2, Fs);

    % reduce the lfp signal to the period of stimulus presentation
    LFP_proc = interp1(t_lfp, filt2, time);
    
    % we are only interested in the stimulus induced fluctuations
%     LFP_proc = LFP_proc - nanmean(LFP_proc(time <= 0));
         
    % the preprocessed lfp and the corresponding time vector
    ex.Trials(n).LFP_prepro = LFP_proc; 
    ex.Trials(n).LFP_prepro_time = time;
end
ex.time = time;

%%
% remove stimulus-driven component (mean subtraction)
% for "internally-generated signals"
stm = [ex.Trials.(ex.exp.e1.type)];
unistm = unique(stm);
lenu = length(unistm);
lens = length(ex.time);
me = cell(1, lenu);
stmc = zeros(1, lenu);
for s = 1:lenu
    me{s} = zeros(1, lens);    
end
for n = 1:N
    ex.Trials(n).LFP_prepro = ex.Trials(n).LFP_prepro(1:lens);
    ex.Trials(n).iLFP_prepro = ex.Trials(n).LFP_prepro;
    stmidx = unistm == ex.Trials(n).(ex.exp.e1.type);
    me{stmidx} = me{stmidx} + ex.Trials(n).LFP_prepro;
    stmc(stmidx) = stmc(stmidx) + 1;
end
for s = 1:lenu
    me{s} = me{s}/stmc(s);
end
% wnd = 0.003;
for n = 1:N
    % mean subtraction for "internally-generated signal"
    stmidx = unistm == ex.Trials(n).(ex.exp.e1.type);
    ex.Trials(n).iLFP_prepro = ex.Trials(n).LFP_prepro - me{stmidx};
    
%     % spike removal
%     ex.Trials(n).LFP_prepro = remove_spk(...
%         ex.Trials(n).LFP_prepro, ex.Trials(n).LFP_prepro_time, ex.Trials(n).Spikes, wnd);
%     ex.Trials(n).iLFP_prepro = remove_spk(...
%         ex.Trials(n).iLFP_prepro, ex.Trials(n).LFP_prepro_time, ex.Trials(n).Spikes, wnd);
%     
%     % mutli taper function
%     [ex.Trials(n).POW, ex.Trials(n).FREQ] = ...
%         mtspectrumc(ex.Trials(n).LFP_prepro(ex.Trials(n).LFP_prepro_time>=0), params);
%     [ex.Trials(n).iPOW, ex.Trials(n).iFREQ] = ...
%         mtspectrumc(ex.Trials(n).iLFP_prepro(ex.Trials(n).LFP_prepro_time>=0), params);
end

function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials
t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
    'UniformOutput', 0); % time of frame onsets
stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);

% also round the average stimulus duration to 2 digits precision
avgstimdur = round(mean(stimdur), 2);
