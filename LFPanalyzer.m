function [para, tuning] = LFPanalyzer(exinfo, pairtype, internal, varargin)
%% comprehensive LFP analysis
% INPUT: exinfo (but assumes a pair of sessions; ex, exinfo(459))
%        pairtype ... 'drug', 'sc', 'ps', 'ps_drug'
%        internal ... 0: no mean subtraction, 1: with, for internally
%        generated LFP
%
% OUTPUT: para, tuning ... structure containing analized data
% 
% written by Katsuhisa (19.03.18)
% +++++++++++++++++++++++++++++++++++++++

% deal with input arguments
if nargin < 2; pairtype = 'drug'; end
if nargin < 3; internal = 0; end
plot_flag = 0;
save_flag = 0;
j = 3;
while j <= nargin - 1
    switch varargin{j-1}
        case 'plot'
            plot_flag = 1;
            j = j + 1;
        case 'save'
            save_flag = 1;
            j = j + 1;
    end
end
    
% stimulus & drugname
para.pairtype = pairtype;
para.stimulus = exinfo.param1;
if exinfo.isRC
    para.stimulus = 'rc';    
    period = {[0.25, 2]};
else
    period = {[0.25, 0.45]};
end
disp(['stimulus: ' para.stimulus])
para.drugname = exinfo.drugname;
para.id = exinfo.id;
para.monkey = exinfo.monkey;

% load lfp data and preprocessing 
% (POW, FREQ, LFP_prepro, LFP_prepro_time)
switch pairtype
    case 'sc'        
        ex = loadCluster(exinfo.fname, 'loadlfp',1);
        [ex0, ex2] = medianSplit(ex, 'sc');
        name1 = 'low sc'; name2 = 'high sc';
        names = {['lowSC_' exinfo.figname], ['highSC_' exinfo.figname]};
    case 'ps'
        ex = loadCluster(exinfo.fname, 'loadlfp',1);
        [ex0, ex2] = medianSplit(ex, 'pupil');
        name1 = 'small ps'; name2 = 'large ps';
        names = {['smallPS_' exinfo.figname], ['largePS_' exinfo.figname]};
    case 'ps_drug'
        ex = loadCluster(exinfo.fname_drug, 'loadlfp',1);
        [ex0, ex2] = medianSplit(ex, 'pupil');
        name1 = 'small ps (drug)'; name2 = 'large ps (drug)';
        names = {['smallPSdrug_' exinfo.figname], ['largePSdrug_' exinfo.figname]};
    otherwise
        ex0 = loadCluster(exinfo.fname, 'loadlfp',1);
        ex2 = loadCluster(exinfo.fname_drug, 'loadlfp',1);
        name1 = 'base'; name2 = para.drugname;
        names = {['base_' exinfo.figname], [para.drugname '_' exinfo.figname]};
end
exs = unitpair_preprocess({ex0, ex2}, period, internal, names);
ex0 = exs{1}; ex2 = exs{2};
% 
% % filter LFP traces into bands
% ex0 = filterLFP(ex0);
% ex2 = filterLFP(ex2);

% LFP properties across stimulus type
tuning = para;
[para.cond(1).lfpstm, tuning.cond(1).results] = LFPbyStm(ex0);
[para.cond(2).lfpstm, tuning.cond(2).results] = LFPbyStm(ex2);

% visualization
if plot_flag==1
    close all;
    visualizer(para, name1, name2, exinfo.figname, save_flag)
end
