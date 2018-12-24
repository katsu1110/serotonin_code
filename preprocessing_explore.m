function ex = preprocessing_explore(ex)
% exploring preprocessing for LFP
% INPUT: ex ... ex-file with LFP data
%              level ... 0, no preprocessing
%                           1, + only baseline correction
%                           2, + bandpass filtering
%                           3, + detrending
%                           4, + line noise removal
%                           5, + removal of trials with outlier
%
% written by Katsuhisa (27.04.18)
% ++++++++++++++++++++++++++++++++

%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% notch filter variables
Fs = 1000;              % sampling frequency
% notchf = [49 51];     % notch filter frequency1
% notchf2 = [99 101];
% notchf3 = [149 151];
% notchord = 2;           % filter order

bpf = [1 100];     % bandpass filter cutoff frequency
bpord = 3;            % filter order

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.2;

% parameters for chronux toolbox
params = define_params;

% TODO: how to handle files with missing LFP entries (such as ma0014,
% 5.44pm)
if any(cellfun(@isempty, {ex.Trials.LFP}))
%     warning('TODO: in frequAnalysis line 44. How to proceed with empty LFP entries.');
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

% %% generate filters
% define notch filter
% [b_notch,a_notch] = butter(notchord, notchf/(Fs/2), 'stop' );
% [b_notch2,a_notch2] = butter(bpord, notchf2/(Fs/2), 'stop');
% [b_notch3,a_notch3] = butter(bpord, notchf3/(Fs/2), 'stop');
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',Fs);
% 
% % define bandpass filter
[b_bpass,a_bpass] = butter(bpord, bpf/(Fs/2), 'bandpass');
% % [b_bpass_notch,a_bpass_notch] = butter(bpord/2, notchf/(Fs/2), 'bandpass');

%% perform functions on each trial lfp

N = length(ex.Trials);
% on = zeros(1, N); % overall noise
for ind = 1:N
   % timing
    t_frame = ex.Trials(ind).Start - ex.Trials(ind).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(ind).LFP_ts - t_frame(1) ; % time rel:stimulus onset

    % raw LFP
    ex.Trials(ind).preproLFP(1).lfp = ex.Trials(ind).LFP;
    
    % we are only interested in the stimulus induced fluctuations
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    time = time - t_frame(1);
    ex.Trials(ind).preproLFP(2).lfp = ex.Trials(ind).preproLFP(1).lfp...
        - nanmean(ex.Trials(ind).preproLFP(1).lfp(time<=0));
    
    %%% bandpass filter
    ex.Trials(ind).preproLFP(3).lfp = filtfilt(b_bpass, a_bpass, ex.Trials(ind).preproLFP(2).lfp);        
    
    % detrending
    ex.Trials(ind).preproLFP(4).lfp = locdetrend(ex.Trials(ind).preproLFP(3).lfp, Fs);
    
    %% notch filter - filters line noise    
    % apply notch filter
%     ex.Trials(ind).preproLFP(5).lfp = filtfilt(b_notch, a_notch, ex.Trials(ind).preproLFP(4).lfp);
%     ex.Trials(ind).preproLFP(5).lfp = filtfilt(b_notch2, a_notch2, ex.Trials(ind).preproLFP(5).lfp);
%     ex.Trials(ind).preproLFP(5).lfp = filtfilt(b_notch3, a_notch3, ex.Trials(ind).preproLFP(5).lfp);
    
%     % remove 50Hz line noise with regression
% %     ex.Trials(ind).preproLFP(5).lfp = rmlinesc(ex.Trials(ind).preproLFP(4).lfp, params,1,0,50);
    ex.Trials(ind).preproLFP(5).lfp = filtfilt(d, ex.Trials(ind).preproLFP(4).lfp);
    
    % reduce the lfp signal to the period of stimulus presentation    
    for v = 1:5
        ex.Trials(ind).preproLFP(v).lfp = interp1( t_lfp, ex.Trials(ind).preproLFP(v).lfp, time );
        % mutli taper function
        [ex.Trials(ind).preproLFP(v).POW, ex.Trials(ind).preproLFP(v).FREQ] = ...
            mtspectrumc(ex.Trials(ind).preproLFP(v).lfp(time>0.35), params);
    end
    
%     on(ind) = mean(ex.Trials(ind).POW);
     
%     % the preprocessed lfp and the corresponding time vector
%     ex.Trials(ind).LFP_prepro = LFP_proc; 
%     ex.Trials(ind).LFP_prepro_time = time;
%     
%     % the preprocessed lfp and the corresponding time vector during
%     % stimulus presentation
%     ex.Trials(ind).LFP_prepro_stm = LFP_proc(time >= 0); 
%     ex.Trials(ind).LFP_prepro_stmtime = time(time >= 0);

end

ex.time = time; % time corresponding to saved lfp signal
ex.time_stm = time(time >= 0);

% % remove trials with too much noise
% ex.Trials(on > 3*std(on)) = [];
end




function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials


t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
    'UniformOutput', 0); % time of frame onsets
stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);


% also round the average stimulus duration to 2 digits precision
avgstimdur = round(mean(stimdur), 2);

end

