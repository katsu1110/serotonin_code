function [stmseq, spikec] = getSc(ex, dt, seqfield, spkfield)
%%
% get stimulus sequence and corresponding spike counts from ex-file
%
% INPUT: ex ... ex-file
%              dt ... bin-duration (sec; 0.001 in default)
%              seqfield ... 'hdx_seq' or 'or_seq'
%              spkfield ... 'Spikes' or 'oSpikes'                   
%
% OUTPUT: stmseq ... stimulus sequence (trials x time-bin)
%                 spikec ... binary (0 or 1) spike counts (trials x time-bin)
%
% EXAMPLE: [stmseq, spikec] = getSc(ex);
%

if nargin < 2; dt = 0.001; end
if nargin < 3; seqfield = 'hdx_seq'; end
if nargin < 4; spkfield = 'Spikes'; end

% use only completed trials
ex.Trials = ex.Trials(abs([ex.Trials.Reward]) > 0);
ntr = length(ex.Trials);

% refresh rate
rf = ex.setup.refreshRate;

% stimulus duration
stmdur = length(ex.Trials(1).(seqfield))/rf;

% time-bin 
ntb = round(stmdur/dt);

% matrix initialization
stmseq = zeros(ntr, ntb); spikec = zeros(ntr, ntb);

% get data
for i = 1:ntr
    % timing of stimulus
    t_strt = ex.Trials(i).Start - ex.Trials(i).TrialStart;
    
    % stm seq & spike count
    ex.Trials(i).(spkfield) = ex.Trials(i).(spkfield) - t_strt(1);
    begin = [0, 1];
    for t = 1:ntb
%         if sum(ex.Trials(i).(spkfield) >= begin(1) & ...
%             ex.Trials(i).(spkfield) <= begin(1) + dt)
%             spikec(i, t) = 1;
%         end
        spikec(i, t) = sum(ex.Trials(i).(spkfield) >= begin(1) & ...
            ex.Trials(i).(spkfield) <= begin(1) + dt);
        if begin(1) > begin(2)/rf
            begin(2) = begin(2) + 1;
        end
        stmseq(i, t) = ex.Trials(i).(seqfield)(begin(2));
        begin(1) = begin(1) + dt;
    end
end

