function si = getSI(fname, analysiswnd, binsize, Fs)
%%
% compute synchronization index 
%
% INPUT: fname ... c1 data (single-unit)
%             analysiswnd ... analysis window in a trial
% OUTPUT: si ... synchronization index (trials x SU, MU)
%
% EXAMPLE: si = getSI('Z:\data\kaki\0256\ka_0256_c1_sortLH_12.32.grating.ORxRC.mat', [0.8 0])
%

if nargin < 2; analysiswnd = [0.8 2]; end
if nargin < 3; binsize = 50; end
if nargin < 4; Fs = 1000; end

% load c1
load(fname, 'ex')
ex1 = ex;

% load c0
load(strrep(fname, 'c1', 'c0'), 'ex')
ex0 = ex;

if ex1.exp.StimPerTrial == 4
    stmdur = 0.45;
elseif ex1.exp.StimPerTrial == 1
    stmdur = 2;
end

% only completed trials
comp = abs([ex1.Trials.Reward]) > 0;
ex1.Trials = ex1.Trials(comp);
ex0.Trials = ex0.Trials(comp);

% spike count
spk1 = getSpks([ex1.Trials], analysiswnd);
spk0 = getSpks([ex0.Trials], analysiswnd);

% synchronization index
ntr = length(ex1.Trials);
si = nan(ntr, 2);
for n = 1:ntr
    si(n, 1) = spktrain2psi(spk1{n}, [analysiswnd(1) stmdur], binsize, Fs);
    si(n, 2) = spktrain2psi(spk0{n}, [analysiswnd(1) stmdur], binsize, Fs);
end