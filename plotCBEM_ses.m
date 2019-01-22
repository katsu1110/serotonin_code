function plotCBEM_ses(CBEM_fit)

%% 
% plot filters 
figure(1);
clf;
subplot(1,3,1);
hold on
plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{1}(1:CBEM_fit.stimNumBasisVectors));
yy1 = get(gca, 'YLim');

ylabel('filter weight');
xlabel('time (ms)');
title('excitatory filter');
hold off


subplot(1,3,2);
hold on
plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{2}(1:CBEM_fit.stimNumBasisVectors));
yy2 = get(gca, 'YLim');
ylabel('filter weight');
xlabel('time (ms)');
title('inhibitory filter');
hold off


subplot(1,3,3);
hold on
plot((1:size(CBEM_fit.spkHistBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.spkHistBasisVectors*CBEM_fit.h_spk(1:CBEM_fit.spkHistNumBasisVectors));
% yy3 = get(gca, 'YLim');
ylabel('filter weight');
xlabel('time (ms)');
title('spike history filter');
hold off

yy = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
for i = 1:2
    subplot(1,3,i)
    ylim(yy)
end