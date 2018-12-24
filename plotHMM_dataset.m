function plotHMM_dataset
% plot results of 'hmms = fitHMM_dataset'

%%
% load data
if mean(ismember('gpfs0', cd))==1
    load('/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/hmms10.mat') 
else
   load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\hmms10.mat')
end

disp(['the number of original sessions: ' num2str(length(hmms.session))])
pass = zeros(1,length(hmms.session));
for i = 1:length(hmms.session)
    if hmms.session(i).estimate_drug(2).exist == 1
        if min([hmms.session(i).estimate(1).fit.likelihood, hmms.session(i).estimate(2).fit.likelihood, ...
                hmms.session(i).estimate_drug(1).fit.likelihood, hmms.session(i).estimate_drug(2).fit.likelihood]) > 0.5
            pass(i) = 1;
        end
    end
end
hmms.session = hmms.session(pass==1);
lenses = length(hmms.session);
disp(['the number of analized sessions: ' num2str(lenses)])

ser = zeros(lenses,1);
for i = 1:lenses
    if ~isempty(strfind(hmms.session(i).fname_drug, '5HT'))
        ser(i) = 1;
    end
end

%%
% plot figures
close all;
on_col = summer(1);
off_col = spring(1);

% firing rate in each state
h = figure(1);
fr = zeros(lenses, 6);
for i = 1:lenses
    fr(i,1) = hmms.session(i).estimate(1).fit.fr(1);
    fr(i,2) = hmms.session(i).estimate_drug(1).fit.fr(1);
    fr(i,3) = hmms.session(i).estimate(2).fit.fr(1);
    fr(i,4) = hmms.session(i).estimate(2).fit.fr(2);
    fr(i,5) = hmms.session(i).estimate_drug(2).fit.fr(1);
    fr(i,6) = hmms.session(i).estimate_drug(2).fit.fr(2);
end
subplot(2,3,1)
unity_scatter(fr(ser==0, 1), fr(ser==0, 2))
xlabel('control')
ylabel('NaCl')
subplot(2,3,2)
unity_scatter(fr(ser==0, 3), fr(ser==0, 4))
xlabel('state 1 (baseline)')
ylabel('state 2 (baseline)')
subplot(2,3,3)
unity_scatter(fr(ser==0, 5), fr(ser==0, 6))
xlabel('state 1 (NaCl)')
ylabel('state 2 (NaCl)')
subplot(2,3,4)
unity_scatter(fr(ser==1, 1), fr(ser==1, 2))
xlabel('control')
ylabel('5HT')
subplot(2,3,5)
unity_scatter(fr(ser==1, 3), fr(ser==1, 4))
xlabel('state 1 (baseline)')
ylabel('state 2 (baseline)')
subplot(2,3,6)
unity_scatter(fr(ser==1, 5), fr(ser==1, 6))
xlabel('state 1 (5HT)')
ylabel('state 2 (5HT)')
set(h, 'Name', 'firing rate in each state', 'NumberTitle', 'off')

% episode duration
h = figure(2);
mesd = zeros(lenses, 8);
for i = 1:lenses
    mesd(i,1) = mean(hmms.session(i).estimate(2).fit.duration.state(1).duration);
    mesd(i,2) = std(hmms.session(i).estimate(2).fit.duration.state(1).duration);
    mesd(i,3) = mean(hmms.session(i).estimate(2).fit.duration.state(2).duration);
    mesd(i,4) = std(hmms.session(i).estimate(2).fit.duration.state(2).duration);
    mesd(i,5) = mean(hmms.session(i).estimate_drug(2).fit.duration.state(1).duration);
    mesd(i,6) = std(hmms.session(i).estimate_drug(2).fit.duration.state(1).duration);
    mesd(i,7) = mean(hmms.session(i).estimate_drug(2).fit.duration.state(2).duration);
    mesd(i,8) = std(hmms.session(i).estimate_drug(2).fit.duration.state(2).duration);
end
subplot(2,2,1)
scatter(mesd(ser==0, 1), mesd(ser==0, 5), 30, 'o', 'markerfacecolor', off_col, ...
    'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
xlabel('baseline')
ylabel('NaCl')
title('state 1')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,2,2)
scatter(mesd(ser==0, 3), mesd(ser==0, 7), 30, 'o', 'markerfacecolor', on_col, ...
    'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
xlabel('baseline')
title('state 2')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,2,3)
scatter(mesd(ser==1, 1), mesd(ser==1, 5), 30, 'o', 'markerfacecolor', off_col, ...
    'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
xlabel('baseline')
ylabel('5HT')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,2,4)
scatter(mesd(ser==1, 3), mesd(ser==1, 7), 30, 'o', 'markerfacecolor', on_col, ...
    'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
xlabel('baseline')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(h, 'Name', 'episode duration', 'NumberTitle', 'off')
% mean vs std
% subplot(2,3,1)
% scatter(mesd(:, 1), mesd(:, 2), 30, 'o', 'markerfacecolor', off_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 1)')
% ylabel('SD of episode duration (state 1)')
% title('baseline')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% subplot(2,3,2)
% scatter(mesd(:, 3), mesd(:, 4), 30, 'o', 'markerfacecolor', on_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 2)')
% ylabel('SD of episode duration (state 2)')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% subplot(2,3,3)
% scatter(mesd(ser==0, 5), mesd(ser==0, 6), 30, 'o', 'markerfacecolor', off_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 1)')
% ylabel('SD of episode duration (state 1)')
% title('NaCl')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% subplot(2,3,4)
% scatter(mesd(ser==0, 7), mesd(ser==0, 8), 30, 'o', 'markerfacecolor', on_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 2)')
% ylabel('SD of episode duration (state 2)')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% subplot(2,3,5)
% scatter(mesd(ser==1, 5), mesd(ser==1, 6), 30, 'o', 'markerfacecolor', off_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 1)')
% ylabel('SD of episode duration (state 1)')
% title('5HT')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% subplot(2,3,6)
% scatter(mesd(ser==1, 7), mesd(ser==1, 8), 30, 'o', 'markerfacecolor', on_col, ...
%     'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)
% xlabel('mean of episode duration (state 2)')
% ylabel('SD of episode duration (state 2)')
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% set(h, 'Name', 'episode duration', 'NumberTitle', 'off')

% variance explained
h = figure(3);
varexp = zeros(lenses, 2);
for i = 1:lenses
    varexp(i,1) = hmms.session(i).estimate(2).fit.variance_explained;
    varexp(i,2) = hmms.session(i).estimate_drug(2).fit.variance_explained;
end
subplot(1,2,1)
unity_scatter(varexp(ser==0,1), varexp(ser==0,2))
xlabel('control')
ylabel('NaCl')
title('variance explained')
hold on;
subplot(1,2,2)
unity_scatter(varexp(ser==1,1), varexp(ser==1,2))
xlabel('control')
ylabel('5HT')


% % likelihood of fitting
% li = zeros(lenses, 2);
% for i = 1:lenses
%     li(i,1) = hmms.session(i).estimate(1).fit.likelihood;
%     li(i,2) = hmms.session(i).estimate_drug(1).fit.likelihood;
% end
% subplot(2,2,3)
% unity_scatter(li(ser==0,2), li(ser==0,2))
% xlabel('control')
% ylabel('NaCl')
% hold on;
% subplot(2,2,4)
% unity_scatter(li(ser==1,2), li(ser==1,2))
% xlabel('control')
% ylabel('5HT')
% title('likelihood')
% set(h, 'Name', 'variance explained & likelihood', 'NumberTitle', 'off')

% on-episode duration vs pupil size


% population average of LFP aligned to the on-off switch


% LFP power in each state
