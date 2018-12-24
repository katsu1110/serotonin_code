function hmms = fitHMM_dataset(ex, dt, seqfield, spkfield, wnd)
%% 
% fit HMM to the spike counts
%
% INPUT: ex ... ex-file
%              dt ... bin-duration (sec; 0.001 in default)
%              seqfield ... 'hdx_seq' or 'or_seq'
%              spkfield ... 'Spikes' or 'oSpikes'                   
%              wnd ... analysis window ([0.8 0] in default)
%
% OUTPUT: hmms .... struct
%
%

%%
% X and Y
[X, Y] = getSc(ex, dt, seqfield, spkfield);
stmdur = dt*size(X, 2);
v = 0:dt:dt*(size(X, 2)-1);
range = v >= wnd(1) & v <= (stmdur - wnd(2));
X = X(:, range); Y = Y(:, range);

%%
% fit GPFA
% try
    hmms.gpfa = fitGPFA(Y);
% catch
%     hmms.gpfa = nan;
% end

% fit HMM
% try
    hmms.hmm = fitHMM(Y, 2, 1);
% catch
%     hmms.hmm = nan;
% end
hmms.stm = X;
