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
% fit GPFA
% try
[X, Y] = getSc(ex, dt(1), seqfield, spkfield);
stmdur = dt(1)*size(X, 2);
v = 0:dt(1):dt(1)*(size(X, 2)-1);
range = v >= wnd(1) & v <= (stmdur - wnd(2));
Xgpfa = X(:, range); Y = Y(:, range); 
hmms.gpfa = fitGPFA(Y);
% catch
%     hmms.gpfa = nan;
% end

% fit HMM
% try
[X, Y] = getSc(ex, dt(2), seqfield, spkfield);
stmdur = dt(2)*size(X, 2);
v = 0:dt(2):dt(2)*(size(X, 2)-1);
range = v >= wnd(1) & v <= (stmdur - wnd(2));
Xhmm = X(:, range); Y = Y(:, range);
Y(Y > 0) = 1;
hmms.hmm = fitHMM(Y, 2, 1);
% catch
%     hmms.hmm = nan;
% end
hmms.stm_gpfa = Xgpfa;
hmms.stm_hmm = Xhmm;
