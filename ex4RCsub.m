function [stmMat, actMat] = ex4RCsub(ex, stmname, fieldname)
% extract relevant information from ex file to perform a RC subspace
% analysis using 'reverse_corr_subspace.m'
% INPUT:
% ex ... ex file
% stmname ... (e.g. 'hdx', 'or'...)
% fieldname ... name of the field in 'ex.Trials' (e.g. 'Spikes',
% 'LFP_prepro')
% 
% written by Katsuhisa (06.04.18)
% ++++++++++++++++++++++++++++++++++++++

% remove instruction or incomplete trials
trls = ex.Trials;
trls = trls(abs([trls.Reward])>0); % Only trls where monkey responds
if isfield(trls, 'instructionTrial')
    trls = trls([trls.instructionTrial]~=1);
end
ntr = length(trls);

% store data into matrix
nframe = length(trls(1).([stmname '_seq']));
stmMat = nan(ntr, nframe);
actMat = zeros(ntr, nframe);
switch fieldname
    case 'Spikes'
        stmtime = linspace(0, ex.fix.duration, nframe);
        for n = 1:ntr
            trls(n).binvec = spk2vec(trls(n).Spikes, stmtime);
        end
    case 'LFP_z'
        stmtime = trls(1).(['LFP_prepro' '_time']) > 0 & trls(1).(['LFP_prepro' '_time']) <= ex.fix.duration;
        for n = 1:ntr
            trls(n).binvec = tcbin(trls(n).(fieldname)(stmtime), nframe);
        end
    otherwise
        stmtime = trls(1).(['LFP_prepro' '_time']) > 0 & trls(1).(['LFP_prepro' '_time']) <= ex.fix.duration;
        for n = 1:ntr
            trls(n).binvec = tcbin(trls(n).(fieldname)(stmtime), nframe);
        end
end
fieldname = 'binvec';
stmtime = 1:nframe;
for n = 1:ntr
    stmMat(n, :) = trls(n).([stmname '_seq']);
    actMat(n, :) = trls(n).(fieldname)(stmtime);
end

% spikes to a vector with a constant length
function spkvec = spk2vec(spk, stmtime)
spkvec = zeros(1, length(stmtime));
for i = 1:length(spk)
    [~, idx] = min(abs(stmtime - spk(i)));
    spkvec(idx) = 1;
end