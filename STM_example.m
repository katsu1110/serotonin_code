function STM_example
%%
% visualize example trace and prediction
%

datapath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data\ka_0278_c1_sortLH_all.grating.ORxRC_5HT\';

prep = load([datapath '\preprocessed_base.mat']);

range = 1000:1200;

close all;
figure;
plot(range, prep.data{2}.calcium(range), '-k')
for i = 1:length(range)
    if prep.data{2}.spikes(range(i)) > 0
        hold on;
        plot(range(i)*[1 1], 0.1*[-1 1] -1.5, '-k')
    end
end
xlim([range(1) range(end)])
set(gca, 'box', 'off', 'tickdir', 'out')
axis off