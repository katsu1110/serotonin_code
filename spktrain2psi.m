function si = spktrain2psi(spkt, analysiswnd, binsize, Fs)
%%
% compute synchronized index (var/mean; fano factor) from a single spike train 
% 

if nargin < 2; analysiswnd = [floor(spkt(1)) ceil(spkt(end))]; end
if nargin < 3; binsize = 50; end
if nargin < 4; Fs = 1000; end

timebin = analysiswnd(1):(binsize/Fs):analysiswnd(end);
lent = length(timebin);
sc = zeros(1, lent-1);
for t = 1:lent-1
    sc(t) = sum(timebin(t) <= spkt & spkt <= timebin(t+1));
end
me = mean(sc); sd = var(sc);
if me==0
    si = nan;
else
    si = sd/me;
end
