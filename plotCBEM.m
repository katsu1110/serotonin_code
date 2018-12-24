function plotCBEM(cbem)
%%
% plot results from 'fitCBEM.m' for serotonin experiment
%
fnames = cbem.list;
lists = [cbem.goodunit', cbem.is5ht', cbem.animal'];

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
cbem.cond(1).results(outs==1) = [];
cbem.cond(2).results(outs==1) = [];
lenses = size(lists, 1);

close all;
figure;
cols = [0 0 0; 1 0 0];
for a = 1:lena
    for k = 1:2
        if a < 3
            cond = find(lists(:,2)==k-1 & lists(:,3)==a-1);
        else
            cond = find(lists(:,2)==k-1);
        end
        for d = 1:2
            mat_e = nan(length(cond), size(cbem.cond(d).results{1}.stimBasisVectors, 1));
            mat_i = nan(length(cond), size(cbem.cond(d).results{1}.stimBasisVectors, 1));
            mat_h = nan(length(cond), size(cbem.cond(d).results{1}.spkHistBasisVectors, 1));
            for c = 1:length(cond)
                mat_e(c, :, :) = cbem.cond(d).results{cond(c)}.stimBasisVectors*cbem.cond(d).results{cond(c)}.k_s{1}...
                    (1:cbem.cond(d).results{cond(c)}.stimNumBasisVectors);
                mat_i(c, :, :) = cbem.cond(d).results{cond(c)}.stimBasisVectors*cbem.cond(d).results{cond(c)}.k_s{2}...
                    (1:cbem.cond(d).results{cond(c)}.stimNumBasisVectors);
                mat_h(c, :, :) = cbem.cond(d).results{cond(c)}.spkHistBasisVectors*cbem.cond(d).results{cond(c)}.h_spk...
                    (1:cbem.cond(d).results{cond(c)}.spkHistNumBasisVectors);
            end
            subplot(lena*2, 3, 1 + 3*(k-1) + 6*(a-1))
            hold on
%             plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{1}(1:CBEM_fit.stimNumBasisVectors));
            me = mean(mat_e, 1);
            sem = std(mat_e, [], 1)/sqrt(size(mat_e, 1));
            fill_between([1:size(mat_e, 2)]*cbem.cond(d).results{1}.dt*1e3, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot([1:size(mat_e, 2)]*cbem.cond(d).results{1}.dt*1e3, me, 'color', cols(d, :))
            
            
            ylabel({animals{a}, drugs{k},'filter weight'});
            if a==1 && k==1 && d==1
                title('excitatory filter');
            end
            hold off


            subplot(lena*2, 3, 2 + 3*(k-1) + 6*(a-1))
            hold on
%             plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{2}(1:CBEM_fit.stimNumBasisVectors));
            me = mean(mat_i, 1);
            sem = std(mat_i, [], 1)/sqrt(size(mat_i, 1));
            fill_between([1:size(mat_i, 2)]*cbem.cond(d).results{1}.dt*1e3, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot([1:size(mat_i, 2)]*cbem.cond(d).results{1}.dt*1e3, me, 'color', cols(d, :))
            xlabel('time (ms)');
            if a==1 && k==1 && d==1
                title('inhibitory filter');
            end
            hold off


            subplot(lena*2, 3, 3 + 3*(k-1) + 6*(a-1))
            hold on
%             plot((1:size(CBEM_fit.spkHistBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.spkHistBasisVectors*CBEM_fit.h_spk(1:CBEM_fit.spkHistNumBasisVectors));
            me = mean(mat_h, 1);
            sem = std(mat_h, [], 1)/sqrt(size(mat_h, 1));
            fill_between([1:size(mat_h, 2)]*cbem.cond(d).results{1}.dt*1e3, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot([1:size(mat_h, 2)]*cbem.cond(d).results{1}.dt*1e3, me, 'color', cols(d, :))

            if a==1 && k==1 && d==1
                title('spike history filter');
            end
            hold off
        end
    end
end