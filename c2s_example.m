function c2s_example
%%
% visualize example trace and prediction
%

datapath = 'Z:\Katsuhisa\serotonin_project\LFP_project\Data\c2s\data\ka_0173_c1_sortLH_all.grating.ORxRC_NaCl\';

prep = load([datapath '\preprocessed_base.mat']);
prep = prep.data;
pred = load([datapath '\predicted_base.mat']);
pred = pred.data;

close all;
figure;
cols = {[.8 0.6 .3], [.3 0.8 .8]};
thre = 0;
movebin = 200;
x = 1:movebin;
start = 1;
while start + movebin < length(prep{1}.calcium)  
    % trial
    range = start:start+movebin-1;

    % original data
    yorig = double(prep{1}.spikes(range));
    % model prediction
    ypred = pred{1}.predictions(range);
    % correlation
    r = corrcoef(yorig, ypred);
    if thre < r(1,2)
        % update threshold
        thre = r(1,2);

        % plot
        clf
        plot(x, prep{1}.calcium(range), '-', 'color', cols{1})
        hold on;
        yorig = yorig/max(yorig);
        plot(x, yorig -2, '-', 'color', cols{1})
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
        plot(x, ypred -2, '-', 'color', cols{2})
        set(gca, 'box', 'off', 'tickdir', 'out')
        xlim([x(1) x(end)])
        axis off
    end
    start = start + movebin;
end
text(-40, 0, 'LFP', 'fontsize', 6, 'color', cols{1})
text(-40, -2, 'original spike rate', 'fontsize', 6, 'color', cols{1})
text(-40, -3, 'predicted spike rate', 'fontsize', 6, 'color', cols{2})
xlim([-40 x(end)])