function [ex_sps, ex_lps, ex] = medianSplit(ex, index)
%%
% perform median split on ex-file by pupil size and spike count
% INPUT: ex ... ex-file
%        index ... 'pupil', 'sc'
%
% OUTPUT: ex_sps ... small pupil size or fr
%         ex_lps ... large pupil size or fr
% 
% written by Katsuhisa (28.09.17)
% +++++++++++++++++++++++++++++++++++++

if nargin < 2; index = 'pupil'; end

% label the trials
label_seq = label4StmPerTr(ex);
isRC = 0;
maxlab = max(label_seq);
if maxlab < 2
    isRC = 1;
end
comptr = find(label_seq > 0);

% old time or new time of ex-file structure
if isfield(ex.Trials(1),'TrialStart_remappedGetSecs')
      time = 'N';   % new
elseif ~isfield(ex.Trials(1),'TrialStart_remappedGetSecs')...
        && isfield(ex.Trials(1),'TrialStartDatapixx')...
        && isfield(ex.Trials(1).times, 'lastEyeReadingDatapixx')
      time = 'O';   % old
end

% trial loop for timing info
len_tr = length(ex.Trials);
for i = 1:len_tr
    % n-th stimulus
    ex.Trials(i).labelseq = label_seq(i);
    ex.Trials(i).(index) = nan;
    ex.Trials(i).stpos = nan;
    ex.Trials(i).enpos = nan;
    if label_seq(i) > 0
        % timing of start and end of stimulus presentation
        if strcmp(time, 'N')            
            t = ex.Trials(i).Eye.t(1:ex.Trials(i).Eye.n) - ex.Trials(i).TrialStartDatapixx;
            st = ex.Trials(i).Start - ex.Trials(i).TrialStart_remappedGetSecs;           
        elseif time == 'O'
            delta = ex.Trials(i).Eye.t(ex.Trials(i).Eye.n) - ex.Trials(i).TrialStartDatapixx - ex.Trials(i).times.lastEyeReading;
            t = ex.Trials(i).Eye.t(1:ex.Trials(i).Eye.n)-ex.Trials(i).TrialStartDatapixx-delta;
            st = ex.Trials(i).Start - ex.Trials(i).TrialStart;
        end

        % get the timing of start and end of stimulus
        [~,ex.Trials(i).stpos] = min(abs(t-st(1)));
        [~,ex.Trials(i).enpos] = min(abs(t-st(end)));     
    end
end

% index
ncomptr = length(comptr);
val = nan(1, len_tr);
switch index
    case 'pupil'        
        % pupil size time-course during the stimulus presentation period
        for i = 1:ncomptr
            temp = nanmean([ex.Trials(comptr(i)).Eye.v(3,:); ex.Trials(comptr(i)).Eye.v(6,:)],1);
            temp = temp(~isnan(temp) & ~isinf(temp));    
            ex.Trials(comptr(i)).(index) = temp;    
            val(comptr(i)) = nanmean(temp(ex.Trials(comptr(i)).stpos:ex.Trials(comptr(i)).enpos));

            % 'trial'
            if ex.Trials(i).labelseq == maxlab && isRC==0 && i > 3
                val(comptr(i-3)) = val(comptr(i));
                val(comptr(i-2)) = val(comptr(i));
                val(comptr(i-1)) = val(comptr(i));
            end
        end
        % filter the pupil size
        ps = [];
        mseq = HiPaFi(val, 5);      % prepare high-pass filter (my function)
        for i = 1:ncomptr
            ex.Trials(comptr(i)).(index) = ex.Trials(comptr(i)).(index) - mseq(i);                    % high-pass filtering
            ex.Trials(comptr(i)).(index) = LoPaFi(ex.Trials(comptr(i)).(index), 1);      % low-pass filtering (my function) 
            ps = [ps, ex.Trials(comptr(i)).(index)(ex.Trials(comptr(i)).stpos:ex.Trials(comptr(i)).enpos)];
        end

        % z-scoring and extract pupil value
        me = nanmean(ps);
        sd = nanstd(ps);
        for i = 1:ncomptr           
            ex.Trials(comptr(i)).(index) = (ex.Trials(comptr(i)).(index) - me)/sd; 
        end
        % pupil metric
        for i = 1:ncomptr  
            if isRC==1 % last 1/8 (250ms) --- RC
                l = length(ex.Trials(comptr(i)).stpos:ex.Trials(comptr(i)).enpos);
                val(comptr(i)) = nanmean(ex.Trials(comptr(i)).(index)...
                    (ex.Trials(comptr(i)).enpos-round(l/8)+1:ex.Trials(comptr(i)).enpos));
    %             [p_min, t_min] = min(ex.Trials(i).pupil_z);
    %             p_max = max(ex.Trials(i).pupil_z(t_min+1:end));
    %             ex.Trials(i).pupil_val = p_max - p_min;
            else % 4 stimuli per trial
                val(comptr(i)) = nanmean(ex.Trials(comptr(i)).(index));
                % 'trial'
                if ex.Trials(i).labelseq == maxlab && isRC==0 && i > 3
                    val(comptr(i-3)) = val(comptr(i));
                    val(comptr(i-2)) = val(comptr(i));
                    val(comptr(i-1)) = val(comptr(i));
                end
            end
        end
    case 'sc'
        for i = 1:ncomptr           
            % spike count
            val(comptr(i)) = sum(ex.Trials(comptr(i)).Spikes > 0 & ex.Trials(comptr(i)).Spikes < 2);
        end
end

% remove trials with nan (non-completed)
ex.Trials = ex.Trials(~isnan(val));
% median split based on the pupil_val on each stimulus type
[stimparam, stmvals] = getStimParam(ex);
sps_tr = [];
lps_tr = [];
for s = 1:length(stmvals)
    idx = find([ex.Trials.(stimparam)] == stmvals(s));
    lenidx = length(idx);
    [~, sorted_idx] = sort(val(idx), 'descend');    
    lps_tr = [lps_tr, idx(sorted_idx(1:floor(lenidx/2)))];
    sps_tr = [sps_tr, idx(sorted_idx(ceil(lenidx/2):end))];
end
ex_sps = ex;
ex_sps.Trials = ex.Trials(sps_tr);
ex_lps = ex;
ex_lps.Trials = ex.Trials(lps_tr);