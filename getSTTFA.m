function sttfa = getSTTFA(lfpspe, lfptime, spktime, wnd, fs)
%%
% generic function to get spike-triggered averaging LFP spectrogram
% 
% INPUT: lfpspe ... spectrogram (freq x time)
%

if isempty(spktime)
    sttfa = nan;
    return
end

% initialization
nspk = length(spktime);
ncol = length(-wnd:(1/fs):wnd);
sttfa = zeros(size(lfpspe, 1), ncol, nspk);
lent = size(lfpspe, 2);

% STA
for i = 1:nspk
    % time range
    tspk = find(lfptime <= spktime(i), 1, 'last');
    tstrt = tspk - wnd*fs;
    tend = tstrt + ncol - 1;
    
    % take all 
    if tstrt < 1
        delta = abs(tstrt) + 1;
        sttfa(:, :, i) = [repmat(lfpspe(:, 1), 1, delta), lfpspe(:, 1:ncol-delta)];
    elseif tend > lent
        delta = tend - lent;
        sttfa(:, :, i) = [lfpspe(:, end-ncol+delta:end), repmat(lfpspe(:, end), 1, delta-1)];
    else
        sttfa(:, :, i) = lfpspe(:, tstrt:tend);
    end
end