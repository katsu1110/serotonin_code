function [nc, ff] = getNC(fname, analysiswnd)
%%
% compute noise correlation from the RC dataset
% assuming that every stimulus type is identical across trials
%
% INPUT: fname ... c1 data (single-unit)
%             analysiswnd ... analysis window in a trial
% OUTPUT: nc ... noise correlation, ff ... fano factor
%
% EXAMPLE: nc = getNC('Z:\data\kaki\0256\ka_0256_c1_sortLH_12.32.grating.ORxRC.mat', [0.2 0])
%

if nargin < 2; analysiswnd = [0.2 0]; end

% load c1
load(fname, 'ex')
ex1 = ex;

% load c0
load(strrep(fname, 'c1', 'c0'), 'ex')
ex0 = ex;

% only completed trials
comp = abs([ex1.Trials.Reward]) > 0;
ex1.Trials = ex1.Trials(comp);
ex0.Trials = ex0.Trials(comp);

% spike count
[~, spkc1] = getSpks([ex1.Trials], analysiswnd);
[~, spkc0] = getSpks([ex0.Trials], analysiswnd);

% noise correlation
rr = corrcoef(zscore(spkc1), zscore(spkc0));
nc = rr(1, 2);

% fano factor
ff = nanvar(spkc1)/nanmean(spkc1);