function mps = MPspectrogram(ex, exn)
% sort spectrogram obtained by MP as a function of stimulus type
%

lenf = length(exn.freq);
lent = length(exn.MPtime);

mps.time = exn.MPtime;
mps.freq = exn.freq;
mps.stm = unique([ex.Trials.(ex.exp.e1.type)]);
lens = length(mps.stm);
mps.signal = zeros(lent, lens);
mps.energy = zeros(lenf, lent, lens);
Trials = exn.Trials;
for s = 1:lens
    str = find([ex.Trials.(ex.exp.e1.type)] == mps.stm(s));
    nstr = length(str);
    for k = 1:nstr
        mps.signal(:, s) = mps.signal(:, s) + Trials(str(k)).signal';
        mps.energy(:, :, s) = mps.energy(:, :, s) + Trials(str(k)).energy;
    end
    mps.signal(:, s) = mps.signal(:, s)/nstr;
    mps.energy(:, :, s) = mps.energy(:, :, s)./nstr;
end
