function [psintr] = pupil_interaction(exinfo, datatype)
% examine any interaction between pupil size and MU, SU and LFP
% INPUT: exinfo ... single pair of experiment, datatype...'MU','SU'
% OUTPUT: psintr ... output structure


%%
addpath(genpath('Z:\Katsuhisa\serotonin_project\'))
addpath(genpath('Z:\Katsuhisa\code\integrated\'))
%%
% structure initialization
psintr.pupil_timecourse_cntr = [];
psintr.pupil_timecourse_drug = [];
psintr.inter_table = nan(2,2);
psintr.b = nan(3,3);
psintr.perf = nan(1,4);
psintr.drugname = exinfo.drugname;

% load data ---------------------------------
switch datatype
    case {'MU','mu'}
        fname = strrep(exinfo.fname, exinfo.cluster, 'c0'); % cluster 0 <=> multiunit activity
        fname_drug = strrep(exinfo.fname_drug, exinfo.cluster, 'c0'); % cluster 0 <=> multiunit activity
        psintr.datatype = 'mu';
    case {'SU','su'}
        fname = exinfo.fname;
        fname_drug = exinfo.fname_drug;
        psintr.datatype = 'su';
    otherwise
        error('The specified datatype is not supported. Try "MU" or "SU".')
end
% control
ex0 = loadCluster(fname, 'ocul', exinfo.ocul, 'loadlfp', false);
% drug
ex2 = loadCluster(fname_drug, 'ocul', exinfo.ocul, 'loadlfp', false);

% remove fields in ex2 which ex0 does not possess
fields = fieldnames(ex2.Trials);
for i = 1:length(fields)
    if ~isfield(ex0.Trials, fields{i})
        ex2.Trials = rmfield(ex2.Trials, fields{i});
    end
end
fields = fieldnames(ex0.Trials);
for i = 1:length(fields)
    if ~isfield(ex2.Trials, fields{i})
        ex0.Trials = rmfield(ex0.Trials, fields{i});
    end
end

% preprocess pupil data ----------------------------
[ex0_sps, ex0_lps, ex0] = pupilSplit(ex0);
[ex2_sps, ex2_lps, ex2] = pupilSplit(ex2);

% trial number
len_tr0 = length(ex0.Trials);
len_tr2 = length(ex2.Trials);    

% store pupil size time-course ---------------------
% -- drug -- label_tr -- TC ---
nmin = 1000;
for i = 1:len_tr0
    if length(ex0.Trials(i).pupil_z) < nmin
        nmin = length(ex0.Trials(i).pupil_z);
    end
end
for i = 1:len_tr2
    if length(ex2.Trials(i).pupil_z) < nmin
        nmin = length(ex2.Trials(i).pupil_z);
    end
end

psmat0 = zeros(len_tr0 , 2 + nmin);
for i = 1:len_tr0
    psmat0(i,2) = ex0.Trials(i).n_stm;
    psmat0(i,3:end) = ex0.Trials(i).pupil_z(1:nmin);
end
psmat2 = ones(len_tr2 , 2 + nmin);
for i = 1:len_tr2
    psmat2(i,2) = ex2.Trials(i).n_stm;
    psmat2(i,3:end) = ex2.Trials(i).pupil_z(1:nmin);
end
psintr.pupil_timecourse_cntr = psmat0;
psintr.pupil_timecourse_drug = psmat2;

% mean firing rate in the baseline condition  -----------------------  
if exinfo.isRC==1
    bfr = mean([ex0.Trials.spkRate]);
end

% use only the first stimulus in the DG experiments
if ~exinfo.isRC
    ex0.Trials = ex0.Trials(psmat0(:,2)==1);
    ex2.Trials = ex2.Trials(psmat2(:,2)==1);
    len_tr0 = sum(psmat0(:,2)==4);
    len_tr2 = sum(psmat2(:,2)==4);
end

% make matrices for GLM ----------------------------------
mat0 = zeros(len_tr0, 5);
mat2 = ones(len_tr2, 5);

% baseline
for i = 1:len_tr0        

    % mean firing rate in this stimulus
    try
        if exinfo.isRC==1
            mat0(i,2) = bfr;
        else
            mat0(i,2) = exinfo.ratemn(exinfo.ratepar(1:end-1)==ex0.Trials(i).(exinfo.param1));
        end
    catch
        continue
    end

    % firing rate in this trial
    mat0(i,1) = ex0.Trials(i).spkRate;

    % drug condition
    mat0(i,3) = 0;

    % pupil size
    mat0(i,4) = ex0.Trials(i).pupil_val;

    % interaction ... mat0(i,5) = 0
end

% drug 
for i = 1:len_tr2        

    % mean firing rate in this stimulus
    try
        if exinfo.isRC==1
            mat2(i,2) = bfr;
        else
            mat2(i,2) = exinfo.ratemn(exinfo.ratepar_drug(1:end-1)==ex2.Trials(i).(exinfo.param1));
        end
    catch
        continue
    end

    % firing rate in this trial
    mat2(i,1) = ex2.Trials(i).spkRate;

    % drug condition
    mat2(i,3) = 1;

    % pupil size
    mat2(i,4) = ex2.Trials(i).pupil_val;

    % interaction
    mat2(i,5) = mat2(i,4);
end

% data matrix with cross validation ------------
rng(198912220);
v0 = 1:len_tr0; v2 = 1:len_tr2;
train_idx0 = datasample(v0, round(len_tr0/2),'Replace',false);
train_idx2 = datasample(v2, round(len_tr2/2),'Replace',false);
test_idx0 = v0(~ismember(v0, train_idx0));
test_idx2 = v2(~ismember(v2, train_idx2));
trainmat = [mat0(train_idx0,:); mat2(train_idx2,:)];
testmat = [mat0(test_idx0,:); mat2(test_idx2,:)];
trainmat = trainmat(~any(isnan(trainmat),2),:); 
testmat = testmat(~any(isnan(testmat),2),:); 
trainmat(:,3:end) = zscore(trainmat(:,3:end));
testmat(:,3:end) = zscore(testmat(:,3:end));

% variance explained by stimulus only ----------------
rr = corrcoef([trainmat(:,2); testmat(:,2)], [trainmat(:,1); testmat(:,1)]);
psintr.perf(1) = rr(1,2)^2;

% model fit (for gain change) -------------------------
varexp = zeros(2,3);
beta = nan(3,6);
for cv = 1:2
    switch cv
        case 1
            train = trainmat;
            test = testmat;
        case 2
            train = testmat;
            test = trainmat;
    end
    for i = 1:3
        % fitting
        x_init = zeros(1, i);
        c = @(x) cost(x, train, i);    
        options = optimset('MaxFunEvals', 10000);
        [b1] = fminsearch(c, x_init, options);
        b_temp = nan(3,1);
        b_temp(1:length(b1),:) = b1;
        beta(:,cv+(i-1)*2) = b_temp;

        % variance explained
        mout = mymodel(b1, test, i);
        rr = corrcoef(mout, test(:,1));
        varexp(cv, i) = rr(1,2)^2;
    end
end
psintr.perf(2:4) = mean(varexp, 1);
for i = 1:3
    psintr.b(i,:) = mean(beta(:, [i, i+3]),2)';
end


% % model fit (for additive change)
% resid = abs(mat(:,1) - mout);
% for i = 1:3
%     % fitting
%     mdl = fitglm(mat(:,3:2+i), resid, ...
%         'Distribution','normal', 'link', 'identity', 'Intercept', false);
%     b2 = mdl.Coefficients.Estimate';
%     b = nan(1, 6);
%     b(4:3+length(b2)) = b2;
%     exinfo.psglm.b = [exinfo.psglm.b; b];
% 
%     % variance explained
%     mout_add = predict(mdl, mat(:,3:2+i));
%     exinfo.psglm.perf(i+4) = 1 - (var(abs(mout + mout_add - mat(:,1)))/var(mat(:,1)));
% end

% interaction table
%%%%%%%%%%%%%%%%%%%%%%%%
% FR %% 5HT %%% base %%%
% S-ps      %%%      %%%
% L-ps      %%%      %%%
%%%%%%%%%%%%%%%%%%%%%%%%

med0 = median(mat0(:,4));
med2 = median(mat2(:,4));
psintr.inter_table(1,1) = mean(mat2(mat2(:,4) < med2,1));
psintr.inter_table(2,1) = mean(mat2(mat2(:,4) > med2,1));
psintr.inter_table(1,2) = mean(mat0(mat0(:,4) < med0,1));
psintr.inter_table(2,2) = mean(mat0(mat0(:,4) > med0,1));

% normalize the table by the grand mean
psintr.inter_table = psintr.inter_table/...
    mean(mean(psintr.inter_table));

% reverse correlation analysis
if exinfo.isRC
    psintr.rc = RC_PS(ex0, ex2, ex0_sps, ex0_lps, ...
        ex2_sps, ex2_lps, {'Spikes'}, psintr.drugname);
end

% model fit for GLM with pupil size =====================
function mout = mymodel(x, mat, incl)
% fr(t) = fr(s)*exp(sum of weighted predictors)
switch incl
    case 1
        mout = mat(:,2).*exp(x(1)*mat(:,3));
    case 2
        mout = mat(:,2).*exp(x(1)*mat(:,3) + x(2)*mat(:,4));
    case 3
        mout = mat(:,2).*exp(x(1)*mat(:,3) + x(2)*mat(:,4) + x(3)*mat(:,5));
end

function f = cost(x, mat, incl)
mout = mymodel(x, mat, incl);
f = 0;
for i = 1:size(mat,1)
    if mout(i) > 0
        f = f - (mat(i,1)*log(mout(i)) - mout(i));
    end
end