function lookinto_coh
% seemingly there are some outliers in coherence analysis...
%

% load data
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\coh_rc.mat')

% outlier detection
lenses = size(coh{1}, 1);
lenf = size(coh{1}, 2);
for i = 1:lenses
   if max(coh{1}(i, round(lenf/3):end, 2)) > 0.45
       disp(i)
   end
end