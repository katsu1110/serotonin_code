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

% data extraction ============================================
mat_e_ses = nan(lenses, size(cbem.cond(1).results{1}.stimBasisVectors, 1), 2);
mat_i_ses = nan(lenses, size(cbem.cond(1).results{1}.stimBasisVectors, 1), 2);
mat_h_ses = nan(lenses, size(cbem.cond(1).results{1}.spkHistBasisVectors, 1), 2);
for i = 1:lenses
    for d = 1:2
        mat_e_ses(i, :, d) = cbem.cond(d).results{i}.stimBasisVectors*cbem.cond(d).results{i}.k_s{1}...
            (1:cbem.cond(d).results{i}.stimNumBasisVectors);
        mat_i_ses(i, :, d) = cbem.cond(d).results{i}.stimBasisVectors*cbem.cond(d).results{i}.k_s{2}...
            (1:cbem.cond(d).results{i}.stimNumBasisVectors);
        mat_h_ses(i, :, d) = cbem.cond(d).results{i}.spkHistBasisVectors*cbem.cond(d).results{i}.h_spk...
            (1:cbem.cond(d).results{i}.spkHistNumBasisVectors);
    end
end

% visualization ===============================================
close all;
cols = [0 0 0; 1 0 0];

f = 1;

% average time-course
figure(f);
for a = 1:lena % animals
    for k = 1:2 % drug types
        if a < 3 
            cond = find(lists(:,2)==k-1 & lists(:,3)==a-1);
        else
            cond = find(lists(:,2)==k-1);
        end
        for d = 1:2 % base or drug
            mat_e = squeeze(mat_e_ses(cond, :, d));
            mat_i = squeeze(mat_i_ses(cond, :, d));
            mat_h = squeeze(mat_h_ses(cond, :, d));
            
            subplot(lena*2, 3, 1 + 3*(k-1) + 6*(a-1))
            hold on
%             plot((1:size(CBEM_fit.stimBasisVectors,1))*CBEM_fit.dt*1e3  ,CBEM_fit.stimBasisVectors*CBEM_fit.k_s{1}(1:CBEM_fit.stimNumBasisVectors));
            me = mean(mat_e, 1);
            sem = std(mat_e, [], 1)/sqrt(size(mat_e, 1));
            xr = [1:size(mat_e, 2)]*cbem.cond(d).results{1}.dt*1e3;
            fill_between(xr, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot(xr, me, 'color', cols(d, :))            
            xlim([xr(1) xr(end)])
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
            xr = [1:size(mat_i, 2)]*cbem.cond(d).results{1}.dt*1e3;
            fill_between(xr, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot(xr, me, 'color', cols(d, :))
            xlim([xr(1) xr(end)])
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
            xr = [1:size(mat_h, 2)]*cbem.cond(d).results{1}.dt*1e3;
            fill_between(xr, me-sem, me+sem, cols(d, :), 0.5)
            hold on;
            plot(xr, me, 'color', cols(d, :))
            xlim([xr(1) xr(end)])
            if a==1 && k==1 && d==1
                title('spike history filter');
            end
            hold off
        end
    end
end
f = f + 1;

% batch
nrow = floor(sqrt(lenses));
ncol = ceil(lenses/nrow);
for i = 1:lenses
    for d = 1:2
        % color
        if d==1
            col = 'k';
        else
            if lists(i, 2) == 0
                col = 'b';
            else
                col = 'r';
            end
        end
        
        % E filter
        figure(f);
        subplot(nrow, ncol, i)
        hold on
        xr = [1:size(mat_e_ses, 2)]*cbem.cond(d).results{i}.dt*1e3;
        plot(xr, squeeze(mat_e_ses(i, :, d)), 'color', col)
        xlim([xr(1) xr(end)])
        
        % I filter
        figure(f+1);
        subplot(nrow, ncol, i)
        hold on        
        xr = [1:size(mat_i_ses, 2)]*cbem.cond(d).results{i}.dt*1e3;
        plot(xr, squeeze(mat_i_ses(i, :, d)), 'color', col)
        xlim([xr(1) xr(end)])
        
        % H filter
        figure(f+2);
        subplot(nrow, ncol, i)
        hold on
        xr = [1:size(mat_h_ses, 2)]*cbem.cond(d).results{i}.dt*1e3;
        plot(xr, squeeze(mat_h_ses(i, :, d)), 'color', col)
        xlim([xr(1) xr(end)])
    end
    for j = 1:3
        figure(f+j-1);
        subplot(nrow, ncol, i)
        title(fname2title(fnames{i}{2}), 'fontsize', 7)
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
end
figure(f);
set(gcf, 'Name', 'excitatory filter', 'NumberTitle', 'off')
figure(f+1);
set(gcf, 'Name', 'inhibitory filter', 'NumberTitle', 'off')
figure(f+2);
set(gcf, 'Name', 'history filter', 'NumberTitle', 'off')
f = f + 3;

function tlab = fname2title(fname)
dotpos = strfind(fname, '.');
tlab = [fname(1:3) fname(dotpos(end-1)+1:dotpos(end)-1)];