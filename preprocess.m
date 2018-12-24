function ex = preprocess(ex)
%%
% preprocess LFP and eye data
%

%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% stimulus per trials
label_seq = label4StmPerTr(ex);

% notch filter variables
Fs = 1000;              % sampling frequency
notchf = [49 51];     % notch filter frequency
notchord = 2;           % filter order

bpf = [3 90];     % bandpass filter cutoff frequency
bpord = [8 10];            % filter order

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.1;

% generate filters ====================================
% define notch filter
[b_notch,a_notch] = butter(notchord, notchf/(Fs/2), 'stop' );
           
% define bandpass filter for LFPs (Nauhaus et al., 2009)
[b_high, a_high] = butter(bpord(1), bpf(1)/(Fs/2), 'high');
[b_low, a_low] = butter(bpord(2), bpf(2)/(Fs/2), 'low');

% define bandpass filter for ps (Urai et al., 2017)
[bps, aps] = butter(2, [0.01 10]/(round(ex.setup.refreshRate)/2), 'bandpass');

% loop to get data ==================================
ex.Trials = ex.Trials(label_seq > 0);
N = length(ex.Trials);

% initialization
lfps = [];
eyeRX = [];
eyeRY = [];
eyeLX = [];
eyeLY = [];
psR = [];
psL = [];
dpsR = [];
dpsL = [];

lfps_timing = ones(N, 2);
eye_timing = ones(N, 4);
for n = 1:N    
    % LFP ================================
    % timing information for LFP
    t_frame = ex.Trials(n).Start - ex.Trials(n).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(n).LFP_ts - t_frame(1) ; % time rel:stimulus onset
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    ex.Trials(n).LFP_prepro_time = time - t_frame(1);
        
    % reduce the lfp signal to the period of stimulus presentation
    lfps_temp = interp1(t_lfp, ex.Trials(n).LFP, time);
    lfps = [lfps, lfps_temp];
    
    % start and end of the trial
    if n > 1
        lfps_timing(n, 1) = lfps_timing(n-1, 2) + 1;
    end
    lfps_timing(n, 2) = lfps_timing(n, 1) + length(lfps_temp) - 1;
    
    % Eye ==============================
    t = ex.Trials(n).Eye.t(1:ex.Trials(n).Eye.n)-ex.Trials(n).TrialStartDatapixx;
    st = ex.Trials(n).Start - ex.Trials(n).TrialStart_remappedGetSecs;            

    % get the timing of start and end of stimulus
    [~,stpos] = min(abs(t-st(1)));
    [~,enpos] = min(abs(t-st(end))); 
    
    % downsample to the stimulus frequency
    nonan = ~isnan(ex.Trials(n).Eye.v(3, :));
    rx = downsample(ex.Trials(n).Eye.v(1, nonan), round(500/ex.setup.refreshRate));
    ry = downsample(ex.Trials(n).Eye.v(2, nonan), round(500/ex.setup.refreshRate));
    lx = downsample(ex.Trials(n).Eye.v(4, nonan), round(500/ex.setup.refreshRate));
    ly = downsample(ex.Trials(n).Eye.v(5, nonan), round(500/ex.setup.refreshRate));
    psr_temp = ex.Trials(n).Eye.v(3, nonan);
    psl_temp = ex.Trials(n).Eye.v(6, nonan);
    psr = downsample(psr_temp, round(500/ex.setup.refreshRate));
    psl = downsample(psl_temp, round(500/ex.setup.refreshRate));
    dpsr = downsample([diff(psr_temp(1:2)) diff(psr_temp)], round(500/ex.setup.refreshRate));
    dpsl = downsample([diff(psl_temp(1:2)) diff(psl_temp)], round(500/ex.setup.refreshRate));
    
    % store vectors
    eyeRX = [eyeRX, rx];
    eyeRY = [eyeRY, ry];
    eyeLX = [eyeLX, lx];
    eyeLY = [eyeLY, ly];
    psR = [psR, psr];
    psL = [psL, psl];
    dpsR = [dpsR, dpsr];
    dpsL = [dpsL, dpsl];
    
    % timings
    lenv = length(rx);
    if n > 1
        eye_timing(n, 1) = eye_timing(n-1, 2) + 1;
    end
    eye_timing(n, 2) = eye_timing(n, 1) + lenv - 1;
    eye_timing(n, 3) = round((lenv*stpos)/ex.Trials(n).Eye.n);
    eye_timing(n, 4) = round((lenv*enpos)/ex.Trials(n).Eye.n);
end

% filtering ==================================
% LFP =======================
lfps = filter(b_notch, a_notch, lfps);
lfps = filter(b_high, a_high, lfps);
lfps = filter(b_low, a_low, lfps);

% Pupil ======================
psR = filter(bps, aps, psR);
psL = filter(bps, aps, psL);
dpsR= filter(bps, aps, dpsR);
dpsL = filter(bps, aps, dpsL);

% back to ex-file ==============================
for n = 1:N
    ex.Trials(n).LFP_prepro = lfps(lfps_timing(n,1):lfps_timing(n,2));
    ex.Trials(n).Eye_prepro.stpos = eye_timing(n, 3);
    ex.Trials(n).Eye_prepro.enpos = eye_timing(n, 4);
    ex.Trials(n).Eye_prepro.RX = eyeRX(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.RY = eyeRY(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.LX = eyeLX(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.LY = eyeLY(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.psR = psR(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.psL = psL(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.dpsR = dpsR(eye_timing(n, 1):eye_timing(n, 2));
    ex.Trials(n).Eye_prepro.dpsL = dpsL(eye_timing(n, 1):eye_timing(n, 2));
end

% remove unnecessary fields
rmf = {'Eye', 'LFP'};
for i = 1:length(rmf)
    ex.Trials = rmfield(ex.Trials, rmf{i});
end

%%
% remove stimulus-driven component (mean subtraction)
% for "internally-generated signals"
stm = [ex.Trials.(ex.exp.e1.type)];
unistm = unique(stm);
lenu = length(unistm);
lens = length(ex.Trials(end).LFP_prepro_time);
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
for n = 1:N
    % mean subtraction for "internally-generated signal"
    stmidx = unistm == ex.Trials(n).(ex.exp.e1.type);
    ex.Trials(n).iLFP_prepro = ex.Trials(n).LFP_prepro - me{stmidx};
end

function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials
t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
    'UniformOutput', 0); % time of frame onsets
stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);

% also round the average stimulus duration to 2 digits precision
avgstimdur = round(mean(stimdur), 2);
