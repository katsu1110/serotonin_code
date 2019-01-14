function plot_HMM(hmms)
%%
% visualize results from hmm analysis
%
% - on-off duration (average)
% - variance explained
% - mean vs SD of on-off duration
% - time-constant fit by exp

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

%
fnames = hmms.list;
lists = [hmms.goodunit', hmms.is5ht', hmms.animal'];

close all;
disp(['5HT: ' num2str(sum(lists(:,2)==1)) ' pairs'])
disp(['NaCl: ' num2str(sum(lists(:,2)==0)) ' pairs'])

animals = {'kaki', 'mango', 'both'};
drugs = {'NaCl', '5HT'};
lena = length(animals);
lenses = size(lists, 1);

% exclude no data sessions
outs = zeros(lenses, 1);
%     outs(6:9) = 1; % mango
%     outs = lists(:,1);
%     outs = 1*(outs < 1);
fnames(outs==1) = [];
lists(outs==1, :) = [];
hmms.cond(1).results(outs==1) = [];
hmms.cond(2).results(outs==1) = [];
lenses = size(lists, 1);

% data extraction ============================================
mean_dur = nan(lenses, 2, 2);
sd_dur = nan(lenses, 2, 2);
variance_explained = nan(lenses, 2);
tau = nan(lenses, 2, 2);
for i = 1:lenses
    % index for state
    idx = [1, 2];
    if hmms.cond(1).results{i}.hmm.state(1).fr > hmms.cond(1).results{i}.hmm.state(2).fr
        idx = [2, 1];
    end
    for d = 1:2
        % mean off duration
        mean_dur(i, 1, d) = nanmean(hmms.cond(d).results{i}.hmm.state(idx(1)).duration);
        % mean on duration
        mean_dur(i, 2, d) = nanmean(hmms.cond(d).results{i}.hmm.state(idx(2)).duration);
        % SD of off duration
        sd_dur(i, 1, d) = nanstd(hmms.cond(d).results{i}.hmm.state(idx(1)).duration);
        % SD of on duration
        sd_dur(i, 2, d) = nanstd(hmms.cond(d).results{i}.hmm.state(idx(2)).duration);
        % variance explained by HMM
        variance_explained(i, d) = hmms.cond(d).results{i}.hmm.variance_explained;
        % time-constant for off duration
        try
            pd = fitdist(hmms.cond(d).results{i}.hmm.state(idx(1)).duration', 'Exponential');
            tau(i, 1, d) = pd(1).mu;
        catch
            tau(i, 1, d) = 0;
        end
        
        % time-constant for on duration
        try 
            pd = fitdist(hmms.cond(d).results{i}.hmm.state(idx(2)).duration', 'Exponential');
            tau(i, 2, d) = pd(1).mu;
        catch
            tau(i, 2, d) = 0;
        end
    end
end

% visualization ===============================================
close all;
% cols = [0 0 0; 1 0 0];

f = 1;
figure(f);

% average time-course
for a = 1:lena % animals
    for k = 1:2 % drug types
        if a < 3 
            cond = find(lists(:,2)==k-1 & lists(:,3)==a-1);
        else
            cond = find(lists(:,2)==k-1);
        end

        % mean duration (off)
        subplot(lena*2, 5, 1 + 5*(k-1) + 10*(a-1))
        unity_scatter(squeeze(mean_dur(cond, 1, 1)), squeeze(mean_dur(cond, 1, 2)))
        if k==2
            xlabel('base')
            ylabel('drug')
            title('off duration')
        end
        
        % mean duration (on)
        subplot(lena*2, 5, 2 + 5*(k-1) + 10*(a-1))
        unity_scatter(squeeze(mean_dur(cond, 2, 1)), squeeze(mean_dur(cond, 2, 2)))
        if k==2
            xlabel('base')
            ylabel('drug')
            title('on duration')
        end
        
        % mean vs SD
        subplot(lena*2, 5, 3 + 5*(k-1) + 10*(a-1))
        me = [squeeze(mean_dur(cond, 1, 1)); squeeze(mean_dur(cond, 1, 2)); ...
            squeeze(mean_dur(cond, 2, 1)); squeeze(mean_dur(cond, 2, 2))];
        sd = [squeeze(sd_dur(cond, 1, 1)); squeeze(sd_dur(cond, 1, 2)); ...
            squeeze(sd_dur(cond, 2, 1)); squeeze(sd_dur(cond, 2, 2))];
        ncond = length(cond);
        lab = [ones(ncond, 1); 2*ones(ncond, 1); 3*ones(ncond, 1); 4*ones(ncond, 1)];
        ls_scatter(me, sd, lab)
        if k==2
            xlabel('mean duration')
            ylabel('SD duration')
            title('off_{base}, off_{drug}, on_{base}, on_{drug}')
        end
        
        % tau (off)
        subplot(lena*2, 5, 4 + 5*(k-1) + 10*(a-1))
        unity_scatter(squeeze(tau(cond, 1, 1)), squeeze(tau(cond, 1, 2)))
        if k==2
            xlabel('base')
            ylabel('drug')
            title('tau (off)')
        end
        
        % tau (on)
        subplot(lena*2, 5, 5 + 5*(k-1) + 10*(a-1))
        unity_scatter(squeeze(tau(cond, 2, 1)), squeeze(tau(cond, 2, 2)))
        if k==2
            xlabel('base')
            ylabel('drug')
            title('tau (on)')
        end
    end
end

