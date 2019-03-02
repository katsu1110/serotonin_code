function c2s_example
%%
% visualize example trace and prediction
%

datapath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data\ka_0173_c1_sortLH_all.grating.ORxRC_NaCl\';

addpath(genpath(['Z:/Katsuhisa/code/integrated/cbrewer']))
cb = cbrewer('div', 'PiYG', 10);
prep = load([datapath '\preprocessed_base.mat']);
prep = prep.data;
pred = load([datapath '\predicted_base.mat']);
pred = pred.data;

close all;
% cols = {[.8 0.6 .3], [.3 0.8 .8]};
cols = {[cb(end, :); cb(1,:)], [cb(end-2, :); cb(3,:)]};
thre = 0;
movebin = 200;
x = 1:movebin;
for i = 1:2
    start = 1;
    figure;
    while start + movebin < length(prep{i}.calcium)  
        % trial
        range = start:start+movebin-1;

        % original data
        yorig = double(prep{i}.spikes(range));
        % model prediction
        ypred = pred{i}.predictions(range);
        % correlation
        r = corrcoef(yorig, ypred);
        if thre < r(1,2)
            % update threshold
            thre = r(1,2);

            % plot
            clf
            plot(x, prep{i}.calcium(range), '-', 'color', cols{i}(1,:))
            hold on;
            yorig = yorig/max(yorig);
            plot(x, yorig -2, '-', 'color', cols{i}(1,:))
    %         for i = 1:length(yorig)
    %             if yorig(i)==1
    %                 plot(x(i)*[1 1], [-1 -0.8], '-', 'color', cols{1})
    %                 hold on;
    %             end
    %         end
    %         for i = 1:length(ypred)
    % %             if ypred(i) > 0.5
    %             if poissrnd(ypred(i)) >= 1
    %                 plot(x(i)*[1 1], [-1.5 -1.3], '-', 'color', cols{2})
    %                 hold on;
    %             end
    %         end
            ypred = ypred/max(ypred);
            plot(x, ypred -3, '-', 'color', cols{i}(end,:))
            set(gca, 'box', 'off', 'tickdir', 'out')
            xlim([x(1) x(end)])
            axis off
        end
        start = start + movebin;
    end
    text(-40, 0, 'LFP', 'fontsize', 6, 'color', cols{i}(1,:))
    text(-40, -2, 'original spike rate', 'fontsize', 6, 'color', cols{i}(1,:))
    text(-40, -3, 'predicted spike rate', 'fontsize', 6, 'color', cols{i}(2,:))
    xlim([-40 x(end)])
end