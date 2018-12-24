function lfpnew = remove_spk(lfp, lfptime, spktime, wnd)
% remove spike from lfp by a linear interpolation
if isempty(spktime)
    lfpnew = lfp;
    return
end
lfpnew = lfp;
nspk = length(spktime);
for n = 1:nspk
    [~, idx1] = min(abs(lfptime - spktime(n) - wnd));
    [~, idx2] = min(abs(lfptime - spktime(n) + wnd));
    if idx1 > idx2
        idx1ori = idx1;
        idx1 = idx2;
        idx2 = idx1ori;
    end
    if idx1 == idx2
        continue
    end
    lfpnew(idx1:idx2) = interp1([idx1, idx2], lfp([idx1, idx2]), idx1:idx2, 'linear');
end