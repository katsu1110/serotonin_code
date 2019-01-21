function stlfp = getSTA(lfptrace, lfptime, spktime, wnd, fs)
%%
% generic function to get spike-triggered averaging LFP
% note that the input 'lfptrace' must be a 0-mean vector
%

if isempty(spktime)
    stlfp = nan;
    return
end

% 0-mean vector
lfptrace = lfptrace - mean(lfptrace(lfptime >= spktime(1) & lfptime <= spktime(end)));

% initialization
nspk = length(spktime);
ncol = length(-wnd:(1/fs):wnd);
stlfp = nan(nspk, ncol);
lent = length(lfptrace);

% STA
for i = 1:nspk
    % time range
    tspk = find(lfptime <= spktime(i), 1, 'last');
    tstrt = tspk - wnd*fs;
    tend = tstrt + ncol - 1;
    
    % take all 
    if tstrt < 1
        continue
%         delta = abs(tstrt) + 1;
%         seg = [lfptrace(1)*ones(1, delta), lfptrace(1:ncol-delta)];
    elseif tend > lent
        continue
%         delta = tend - lent;
%         seg = [lfptrace(end-ncol+delta:end), lfptrace(end)*ones(1, delta-1)];
    else
        seg = lfptrace(tstrt:tend);
    end
    stlfp(i, :) = seg;
end
stlfp(any(isnan(stlfp), 2), :) = [];