function MP_single(ex, internal, filename)
%%
% MP decomposition preprocessing for a single ex-file with preprocessed LFP data
% INPUT: 
%

% preset parameters
Fs = 1000; % sampling frequency
folderName = '../../../LFP_project/Data/MPprepro/';
max_iter = 100;
wrap = 1; % if 1, it does not work...
L = 4096; % 2^N
if (2/ex.exp.StimPerTrial + 0.2)*Fs < 1024
    L = 1024;
elseif (2/ex.exp.StimPerTrial + 0.2)*Fs < 2048
    L = 2048;
end
signalRange = [1 L];
atomList = [];
T = length(ex.Trials(end).LFP_prepro_time);
ntr = length(ex.Trials);
if internal==0
    lfpfield = 'LFP_prepro';
else
    lfpfield = 'iLFP_prepro';
end
% rmdir([folderName filename]) 

% MP decomposition
freq = 0:Fs/L:100;
lenf = length(freq);
inputSignal = nan(L, ntr);
idxs = ones(1, ntr);
for k = 1:ntr
    % padding and mean subtraction for signal processing
    [z, idxs(k)] = padding(ex.Trials(k).(lfpfield), L);
    inputSignal(:, k) = z - mean(z);
end

% MP =========================
filename = strrep(filename, '\', '/');

% perform Gabor decomposition
importData(inputSignal, folderName, filename, signalRange, Fs);
runGabor(folderName, filename, L, max_iter);

% signal reconstruction
exn.MPtime = ex.Trials(end).LFP_prepro_time;
exn.freq = freq;
MPtrials = struct('signal', [], 'energy', []);
parfor k = 1:ntr
    % load data
    gaborInfo = getGaborData(folderName, filename, 1);

    % signal
    mp_signal = reconstructSignalFromAtomsMPP(...
        gaborInfo{k}.gaborData, L, wrap, atomList);
    mp_signal = mp_signal(idxs(k):idxs(k)+T-1);
    MPtrials(k).signal = mp_signal;

    % energy
    rEnergy = reconstructEnergyFromAtomsMPP(...
        gaborInfo{k}.gaborData, L, wrap, atomList);
    rEnergy = rEnergy(1:lenf, idxs(k):idxs(k)+T-1);
    MPtrials(k).energy = rEnergy;
end
exn.Trials = MPtrials;

% autosave =======================
if internal==0
    save([folderName 'Trials/' filename], 'exn', '-v7.3')
else
    save([folderName 'iTrials/' filename], 'exn', '-v7.3')
end
clc
disp([filename ' saved after MP!'])

% subfunction
function [b, pre] = padding(a, L)
% add padding for signal processing
lena = length(a);
b = zeros(1, L);
pre = floor((L - lena)/2) + 1;
b(1:pre-1) = a(1);
b(pre:pre+lena-1) = a;
b(pre+lena:end) = a(end);