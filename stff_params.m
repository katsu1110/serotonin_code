function [nsc, nov, nff] = stff_params(L, overlap)
nsc = floor(L/18);
nov = floor(nsc/(100/overlap));
% nff = max(256, 2^nextpow2(nsc));
nff = 1024;