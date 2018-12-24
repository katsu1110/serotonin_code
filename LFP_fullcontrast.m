function LFP_fullcontrast(dat)
% extract average lfp signal for responses to 100% contrast stimuli

if nargin>0

    dat = dat(~[dat.isc2]);
    dat(1).lfp = []; dat(1).lfpsem = []; 
    dat(1).lfp_drug = []; dat(1).lfpsem_drug = [];
    dat(1).lfp_diff = [];
    dat(1).lfp_ts = [];
    dat(1).pow = []; dat(1).powsem = [];
    dat(1).pow_drug = []; dat(1).powsem_drug = [];
    dat(1).pow_diff = [];
    dat(1).freq = [];
    dat(1).sta = []; dat(1).sta_drug = [];
    dat(1).sta_sd = []; dat(1).sta_sd_drug = [];
    dat(1).nspk = [];
    dat(1).staco = [];
    dat(1).nspkco = [];
    dat(1).staco_drug = [];
    dat(1).nspkco_drug = [];
    
    parfor i = 1:length(dat)
        dat2(i) = setLFP(dat(i));
    end
    dat=dat2;
    save('Data\rcdat.mat', 'dat');
else
    load('Data\rcdat');
end


dat = dat([dat.id]~=237 & [dat.id]~=263 & [dat.id]~=272 & [dat.id]~=280 & [dat.id]~=303 & [dat.id]~=306 & [dat.id]~=313 & ...
    [dat.id]~=324 & [dat.id]~=313 & [dat.id]~=323  & [dat.id]~=331 & ...
    [dat.id]~=4.5 & [dat.id]~=24.5 & [dat.id]~=125.5 & [dat.id]~=127.5 &[dat.id]~=126.5 & [dat.id]~=137.5 ...
    & [dat.id]~=162.5 & [dat.id]~=172.5);   
dat = dat([dat.is5HT]);
% figure;
% evaluateLFP(dat([dat.gslope]<0.9)); 


for i = 1:length(dat)
    h = figure('Name', dat(i).figname, 'Position', [332 394 1119 589]);
    evaluateSingleUnitLFP(dat(i)); 
    savefig(h, [fullfile('FiguresRC', dat(i).figname), '.fig']);
    close(h);
end


end

%%
function dat = setLFP(dat)

[lfpbase, lfpdrug] = LoadLFPfile( dat );


% spike triggered lfp
exspk = loadCluster(dat.fname, 'ocul', dat.ocul);
[dat.sta, dat.sta_sd, dat.nspk] = spktriglfp( exspk, lfpbase);

exspk_drug = loadCluster(dat.fname_drug, 'ocul', dat.ocul);
[dat.sta_drug, dat.sta_sd_drug, dat.nspk_drug] = spktriglfp( exspk_drug, lfpdrug);


% contrast values
vals = [1001 0.125 0.25 0.5 1]; 


expskco = exspk; lfpco = lfpbase;
expskco_drug = exspk_drug; lfpco_drug = lfpdrug;
% loop through the data set and average for each contrast
for vals_i = 1:length(vals)
    
    idx = 1:length(lfpbase.Trials); %[lfpbase.Trials.(dat.param1)] == vals(vals_i);
    lfp(vals_i, :) = mean(vertcat(lfpbase.Trials(idx).LFP_prepro));
    lfpsem(vals_i, :) = std(vertcat(lfpbase.Trials(idx).LFP_prepro))/sqrt(sum(idx));
    pow(vals_i, :) = median(horzcat(lfpbase.Trials(idx).POW), 2);
    powsem(vals_i, :) = mean(horzcat(lfpbase.Trials(idx).POW), 2)/sqrt(sum(idx));
    
    expskco.Trials = exspk.Trials(idx);
    lfpco.Trials = lfpbase.Trials(idx);
    [staco(vals_i, :), ~, nspkco(vals_i)] = spktriglfp( expskco, lfpco);

    
    idx = 1:length(lfpdrug.Trials); %[lfpdrug.Trials.(dat.param1)] == vals(vals_i);
    lfp_drug(vals_i, :) = mean(vertcat(lfpdrug.Trials(idx).LFP_prepro));
    lfpsem_drug(vals_i, :) = std(vertcat(lfpdrug.Trials(idx).LFP_prepro))/sqrt(sum(idx));
    pow_drug(vals_i, :) = median(horzcat(lfpdrug.Trials(idx).POW), 2);
    powsem_drug(vals_i, :) = mean(horzcat(lfpdrug.Trials(idx).POW), 2)/sqrt(sum(idx));
    
    expskco_drug.Trials = exspk_drug.Trials(idx);
    lfpco_drug.Trials = lfpdrug.Trials(idx);
    [staco_drug(vals_i,:), ~, nspkco_drug(vals_i)] = spktriglfp( expskco_drug, lfpco_drug );

    
end

% normaliize the lfp signal
maxlfp = 1;%max(max(abs(lfp)));
dat.lfp = lfp./maxlfp;
dat.lfpsem = lfpsem;
dat.lfp_drug = lfp_drug./maxlfp;
dat.lfpsem_drug = lfpsem_drug;
dat.lfp_diff = dat.lfp - dat.lfp_drug;
dat.lfp_ts = lfpbase.time;

dat.staco = staco;
dat.nspkco = nspkco;
dat.staco_drug = staco_drug;
dat.nspkco_drug = nspkco_drug;

dat.pow = pow;
dat.pow_drug = pow_drug;
dat.pow_diff = dat.pow_drug./dat.pow;
dat.freq = lfpdrug.Trials(1).FREQ;

end


%%
function evaluateLFP(dat)
%% evaluate LFP by plotting mean and standard error

t = dat(1).lfp_ts; % time vector
t2 = [t fliplr(t)]; % for sem patch
n = length(dat);

vals = [1001 0.125 0.25 0.5 1];

col_base = repmat([0.8 0.6 0.4 0.2 0]', 1,3);
col_5ht = col_base;
col_5ht(:,1) = 1;

%%%% LFP vs time
subplot(2,2,[1 2])
for v_i = 1:length(vals)
    
    fun = @(x) x(v_i,:);
    A(v_i,:) = cellfun(fun, {dat.lfp}, 'UniformOutput', 0);
    B(v_i,:) = cellfun(fun, {dat.lfp_drug}, 'UniformOutput', 0);
    
    %%% Average signal
    lfp_base_avg = mean(vertcat(A{v_i,:}));    % baseline    
    lfp_drug_avg = mean(vertcat(B{v_i,:})); % 5HT
    
    %%% plot results
    % Baseline
    pbase(v_i)=plot(t, lfp_base_avg, 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
    % 5HT
    pdrug(v_i)=plot(t, lfp_drug_avg, 'Color', col_5ht(v_i, :), 'LineWidth', 2,...
        'DisplayName', ['co=' num2str(vals(v_i)) ' 5HT']); hold on;
end


%%% plot deviation (SEM)
for v_i = 1:length(vals)
    % baseline
    lfp_base_sem = std(vertcat(A{v_i,:}), 0, 1)/ sqrt(n);
    fill(t2, [lfp_base_avg+lfp_base_sem, fliplr(lfp_base_avg-lfp_base_sem)], ...
        col_base(v_i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5,...
        'DisplayName', ' ');

    
    % 5HT
    lfp_drug_avg = mean(vertcat(B{v_i,:}));
    lfp_drug_sem = std(vertcat(B{v_i,:}), 0, 1) / sqrt(n);
    fill(t2, [lfp_drug_avg+lfp_drug_sem, fliplr(lfp_drug_avg-lfp_drug_sem)],....
        col_5ht(v_i,:),'EdgeColor', 'none','FaceAlpha', 0.5,...
        'DisplayName', ' ');
end

xlim([-0.1 0.45]);
legend([pbase, pdrug], 'Location', 'bestoutside')

ylabel('filtered LFP (\muV)');
xlabel('time relative to stimulus onset (s)');
title('5HT induced changes in the LFP signal');
crossl




%%%% LFP power
f = dat(1).freq';
f2 = [f fliplr(f)]; % for sem patch

subplot(2,2,3)
for v_i = 1:length(vals)
    
    fun = @(x) x(v_i,:);
    A(v_i,:) = cellfun(fun, {dat.pow}, 'UniformOutput', 0);
    B(v_i,:) = cellfun(fun, {dat.pow_drug}, 'UniformOutput', 0);
    
    %%% Average signal
    pow_base_avg = median(vertcat(A{v_i,:}));    % baseline
    pow_drug_avg = median(vertcat(B{v_i,:})); % 5HT
    
    %%% plot results
    % Baseline
    plot(f, pow_base_avg, 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
    % 5HT
    plot(f, pow_drug_avg, 'Color', col_5ht(v_i, :), 'LineWidth', 2,...
        'DisplayName', ['co=' num2str(vals(v_i)) ' 5HT']); hold on;
end

legend('show', 'Location', 'bestoutside')

%%% plot deviation (SEM)
% for v_i = 1:length(vals)
%     % baseline
%     pow_base_avg = mean(vertcat(A{v_i,:}));
%     pow_base_sem = std(vertcat(A{v_i,:}))/ sqrt(n);
%     fill(f2, [pow_base_avg+pow_base_sem, fliplr(pow_base_avg-pow_base_sem)], ...
%         col_base(v_i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5,...
%         'DisplayName', ' ');
%     
%     % 5HT
%     pow_drug_avg = mean(vertcat(B{v_i,:}));
%     pow_drug_sem = std(vertcat(B{v_i,:})) / sqrt(n);
%     fill(f2, [pow_drug_avg+pow_drug_sem, fliplr(pow_drug_avg-pow_drug_sem)],....
%         col_5ht(v_i,:),'EdgeColor', 'none','FaceAlpha', 0.5,...
%         'DisplayName', ' ');
% end


ylabel('Power (\muV^2)');
xlabel('Frequency (Hz)');
xlim([0 90]); 
set(gca, 'YScale', 'log')


%%% difference in power
subplot(2,2,4)

for v_i = 1:length(vals)
    
    fun = @(x) x(v_i,:);
    A(v_i,:) = cellfun(fun, {dat.pow}, 'UniformOutput', 0);
    B(v_i,:) = cellfun(fun, {dat.pow_drug}, 'UniformOutput', 0);
    
    %%% Average signal
    pow_base_avg = mean(vertcat(A{v_i,:}));    % baseline
    pow_drug_avg = mean(vertcat(B{v_i,:})); % 5HT
    
    %%% plot results
    % Baseline
    plot(f, (pow_base_avg-pow_drug_avg)./(pow_base_avg+pow_drug_avg), 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
end

crossl
ylabel('Modulation index of stimulus induced Power (a-b/a+b) (\muV^2)');
xlabel('Frequency (Hz)');
xlim([10 90]); 

end


%%
function evaluateSingleUnitLFP(dat)
%% evaluate LFP by plotting mean and standard error

t = dat.lfp_ts; % time vector
t2 = [t fliplr(t)]; % for sem patch
n = length(dat);

vals = [1001 0.125 0.25 0.5 1];

col_base = repmat([0.8 0.6 0.4 0.2 0]', 1,3);
col_5ht = col_base;
col_5ht(:,1) = 1;

%%%% ================================== LFP vs time
subplot(2,3,[1 2])
for v_i = 1:length(vals)
    
    %%% plot results
    % Baseline
    pbase(v_i) = plot(t, dat.lfp(v_i,:), 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
    % 5HT
    pdrug(v_i) = plot(t, dat.lfp_drug(v_i,:), 'Color', col_5ht(v_i, :), 'LineWidth', 2,...
        'DisplayName', ['co=' num2str(vals(v_i)) ' 5HT']); hold on;
end


% plot deviation (SEM)
for v_i = 1:length(vals)
    % baseline
    fill(t2, [dat.lfp(v_i,:)+dat.lfpsem(v_i,:)./5, fliplr(dat.lfp(v_i,:)-dat.lfpsem(v_i,:)./5)], ...
        col_base(v_i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.1,...
        'DisplayName', ' ');
    
    % 5HT
    fill(t2, [dat.lfp_drug(v_i,:)+dat.lfpsem_drug(v_i,:)./5, ...
        fliplr(dat.lfp_drug(v_i,:)-dat.lfpsem_drug(v_i,:)./5)],....
        col_5ht(v_i,:),'EdgeColor', 'none','FaceAlpha', 0.1,...
        'DisplayName', ' ');
end

xlim([t(1) t(end)]);
ylabel('averaged LFP +/- 1/5th SEM (\muV)');
xlabel('time relative to stimulus onset (s)');
title('5HT induced changes in the LFP signal');
crossl

legend([pbase pdrug], 'Location', 'bestoutside')


%%%% ================================== LFP power
f = dat(1).freq';
kernel = ones(1,3)/3;
subplot(2,3,4)
for v_i = 1:length(vals)
    
    %%% plot results
    % Baseline
    plot(f, conv(dat.pow(v_i,:), kernel, 'same'), 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
    % 5HT
    plot(f, conv(dat.pow_drug(v_i,:), kernel, 'same'), 'Color', col_5ht(v_i, :), 'LineWidth', 2,...
        'DisplayName', ['co=' num2str(vals(v_i)) ' 5HT']); hold on;
end


ylabel('Power (\muV^2)');
xlabel('Frequency (Hz)');
xlim([0 90]); 
set(gca, 'YScale', 'log')


%%%% ================================== difference in power
subplot(2,3,5);
kernel = ones(1,5)/5;
for v_i = 1:length(vals)
    
    %%% plot differene in power
    plot(f, conv(log(dat.pow(v_i,:))-log(dat.pow_drug(v_i,:)), kernel, 'same'),...
        'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
end

crossl
ylabel('delta Power (base-drug) (\muV^2)');
xlabel('Frequency (Hz)');
xlim([10 90]); 


%%%% ============================== spike triggered lfp - contrast dependent
subplot(2,3,3);
t = -0.40:0.001:0.40;

for v_i = 1:length(vals)

    %%% 5HT
    plot(t, dat.staco_drug(v_i,:), 'Color', col_5ht(v_i, :), 'LineWidth', 2, ...
    'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;

    %%% baseline
    plot(t, dat.staco(v_i,:), 'Color', col_base(v_i, :), 'LineWidth', 2, ...
        'DisplayName', ['co=' num2str(vals(v_i)) ' base']); hold on;
    
    
end


title(sprintf(['stimulus dependent STA \n #spks baseline %1.0f %1.0f %1.0f %1.0f %1.0f \n' ...
    '#spks 5HT: %1.0f %1.0f %1.0f %1.0f %1.0f'], dat.nspkco, dat.nspkco_drug))
ylabel('LFP (\muV)')
crossl;
xlim([t(1) t(end)]);

%%%% ================================== spike triggered average (STA)
subplot(2,3,6);
plot(t, dat.sta, 'Color', 'k'); hold on;
plot(t, dat.sta_drug, 'Color', 'r');
fill([t fliplr(t)], [dat.sta+ dat.sta_sd fliplr(dat.sta-dat.sta_sd)], ...
    'k' ,'EdgeColor', 'none', 'FaceAlpha', 0.5,...
    'DisplayName', ' ');
fill([t fliplr(t)], [dat.sta_drug+ dat.sta_sd_drug fliplr(dat.sta_drug-dat.sta_sd_drug)], ...
    'r' ,'EdgeColor', 'none','FaceAlpha', 0.5,...
    'DisplayName', ' ');

xlabel('time rel to spike (s)')
ylabel('LFP (\muV) +/s SEM')
title(sprintf('baseline (k): %1.0f spikes \n drug (r): %1.0f spikes', dat.nspk, dat.nspk_drug)); 
xlim([t(1) t(end)]);
crossl
end




