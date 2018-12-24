function CBEM_fit = fitCBEM(ex, dt, seqfield, spkfield, figon)
%% 
% fit CBEM (Letimer et al., 2018) to the spike counts
%
% INPUT: ex ... ex-file
%              dt ... bin-duration (sec; 0.001 in default)
%              seqfield ... 'hdx_seq' or 'or_seq'
%              spkfield ... 'Spikes' or 'oSpikes'                   
%
% OUTPUT: stmseq ... stimulus sequence (trials x time-bin)
%                 spikec ... binary (0 or 1) spike counts (trials x time-bin)
%
% EXAMPLE: [stmseq, spikec] = getSc(ex);
%

%%
% X and Y
[X, Y] = getSc(ex, dt, seqfield, spkfield);
X = X'; Y = Y';
X = X(:); Y = Y(:);

%%
% CBEM parameters
nLinearRFs = 10;     %number of stimulus basis functions
numShortSHFilts = 5; %each is 0.4ms long - used for refractory period
nSHfilters = 7+numShortSHFilts;%7+

nInh = 1;
nExc = 1;

useLeakCondExc = false;
usePrior = true;

CBEM = setupCBEMsimple(nLinearRFs,nSHfilters,dt,nInh,nExc,numShortSHFilts,useLeakCondExc,usePrior);

%% 
% Convolves stimulus & spike history with filter basis functions
TT = length(X);
SpikeStim = conv2(X,CBEM.stimBasisVectors);
SpikeStim = SpikeStim(1:TT,:);
spkHist = conv2(Y,CBEM.spkHistBasisVectors);
spkHist = [zeros(1,size(spkHist,2)); spkHist(1:TT-1,:)];

%%
% fitting parameters
%this script assumes a CBEM with 2 stimulus conductances (exc and inh)
%and 1 or 2 leak condutances for 

doFitForGL1    = false; %if true, fits the leak conductance values. Otherwise, keeps them fixed

initV  = -60; %-62

priorWeight_e1 = 1;%5^2;
priorWeight_e2 = 1;%5^2;
priorWeight_i = priorWeight_e2/5;%5^2;
CBEM.prior.comp.m = 0;
CBEM.prior.comp.d = 0;
CBEM.prior.k_s_sig{1} = diag([priorWeight_e1*ones(10,1);0.0]);
CBEM.prior.k_s_sig{2} = diag([priorWeight_i*ones(10,1);0.0]);

CBEM.g_l = 25;
CBEM.log_g_l = log(CBEM.g_l);
totalG = 300 - exp(CBEM.log_g_l);
if(numel(CBEM.E_s) > 3)
    %if using 2 leak conductances (to fit both E_l and g_l) - this is the
    %default setup
    initGl = [(CBEM.E_s(3) - initV) (CBEM.E_s(4) - initV); 1 1]\[-exp(CBEM.log_g_l)*(CBEM.E_l - initV);totalG];
    CBEM.k_s{3} = log(initGl(1));
    CBEM.k_s{4} = log(initGl(2));
    
    CBEM.prior.k_s_sig{3} = 1/10;
    CBEM.prior.k_s_sig{4} = 1/10;
else
    %if using 1 leak conductance (to fit g_l with fixed E_l)
    CBEM.E_l      = initV;
    CBEM.E_s(3) = initV;
    CBEM.k_s{3} = log(totalG);
    
    CBEM.prior.k_s_sig{3} = 0/10;
end

fprintf('Fitting CBEM with single linear conductance...\n');
[CBEM_lin] = fitCBEMwithLinearTransferFR(CBEM,SpikeStim,spkHist,Y);
[~,~,CBEMlin_nll]      = fitCBEMwithLinearTransferFR(CBEM_lin,SpikeStim,spkHist,Y,true);

addOnesColumnToStim = true; %if true, function adds on a column of 1's to the end of SpikeStim


%%
% fitting
CBEM_init = CBEM_lin;
% CBEM_init.k_s{1}(:) = 0;
% CBEM_init.k_s{2}(:) = 0;
CBEM_init.k_s{1}(1:10) =  CBEM_init.k_s{1}(1:10) + randn(10,1)*0;
CBEM_init.k_s{2}(1:10) = -CBEM_init.k_s{1}(1:10) + randn(10,1)*0;
CBEM_init.k_s{1}(end) = 20;
CBEM_init.k_s{2}(end) = 20;
CBEM_init.prior.k_s_sig{1} = diag([priorWeight_e2*ones(10,1);0.0]);
CBEM_init.prior.k_s_sig{2} = diag([priorWeight_i*ones(10,1);0.0]);

fprintf('Fitting CBEM with LN excitatory and inhibitory conductances...\n');
[CBEM_fit] = fitCBEMfull(CBEM_init, SpikeStim,spkHist,Y,addOnesColumnToStim,doFitForGL1);

% [~,CBEMnll]      = fitCBEMfull(CBEM_fit,SpikeStim,     spkHist,     Y,     addOnesColumnToStim,false,true);

fprintf('done.\n');

%% 
% plot filters 
if figon
    figure(1);
    clf;
    subplot(1,3,1);
    hold on
%     plot((1:size(CBEM_true.stimBasisVectors,1))*CBEM_true.dt*1e3,CBEM_true.stimBasisVectors*CBEM_true.k_s{1}(1:CBEM_true.stimNumBasisVectors));
    plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{1}(1:CBEM_fit.stimNumBasisVectors));
    %CBEM_true.k_s{1}(1:10) holds the filter weights, CBEM_true.k_s{1}(11) is a baseline weight 

    ylabel('filter weight');
    xlabel('time (ms)');
    title('excitatory filter');
    legend({'true filter','estimated filter'});
    hold off


    subplot(1,3,2);
    hold on
%     plot((1:size(CBEM_true.stimBasisVectors,1))*CBEM_true.dt*1e3,CBEM_true.stimBasisVectors*CBEM_true.k_s{2}(1:CBEM_true.stimNumBasisVectors));
    plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{2}(1:CBEM_fit.stimNumBasisVectors));

    ylabel('filter weight');
    xlabel('time (ms)');
    title('inhibitory filter');
    legend({'true filter','estimated filter'});
    hold off


    subplot(1,3,3);
    hold on
%     plot((1:size(CBEM_true.spkHistBasisVectors,1))*CBEM_true.dt*1e3,CBEM_true.spkHistBasisVectors*CBEM_true.h_spk(1:CBEM_true.spkHistNumBasisVectors));
    plot((1:size(CBEM_fit.spkHistBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.spkHistBasisVectors*CBEM_fit.h_spk(1:CBEM_fit.spkHistNumBasisVectors));

    ylabel('filter weight');
    xlabel('time (ms)');
    title('spike history filter');
    legend({'true filter','estimated filter'});
    hold off
end

%% Simulating from model fit
[CBEM_fit.sim.Mtsp,CBEM_fit.sim.spks,CBEM_fit.sim.V_fit,CBEM_fit.sim.g_s_fit, CBEM_fit.sim.l_s_fit] = ...
    simulateCBEM(CBEM_fit,SpikeStim,10);