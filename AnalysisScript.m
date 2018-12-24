% compare the static and drifting gratings as well as the depth dependency

clc
clear
close all
%% 1. load all data in the workspace
fdir = 'Z:\data\kaki\0290';
listing = dir(fdir);
listing = listing(3:end); % remove . and ..

% load all ex files.
for i = 1:length(listing)
    ex_all{i} = loadCluster(fullfile(fdir, listing(i).name));
    ex_all{i} = frequAnalysis(ex_all{i}, 'notchf', [28 35], 'notchord', 1);
end

clearvars listing

%% 2. look at one of the raw data

fs = 1000; %sampling frequency
nfft = 1024;

for i = 1:length(ex_all)
    h = figure('Name', ex_all{i}.Header.fileName, ...
        'Position', [296 320 1472 658]);
    cwt( ex_all{i}.LFP.v, 1000);
    savefig(h, [ex_all{i}.Header.fileName(15:end) '.fig'], 'compact');
    delete(h);
end


%%

lfp = ex_all{2}.LFP.v;
[pow, f] = pmtm( lfp, 2, nfft, fs);

figure;
subplot(2,1,1); plot(lfp, 'DisplayName', 'unfiltered'); 
ylabel('\muV'); xlabel('time (unit unknown)');
subplot(2,1,2); plot(f, pow, 'DisplayName', 'unfiltered'); 
xlim([2 100]); ylabel('Power (\muV^2)'); xlabel('Frequency');
% set(gca, 'YScale', 'log');

% filter the signal
[b_notch,a_notch] = butter(4, [10 35]./(fs/2), 'stop' );
[b_bp,a_bp] = butter(2, [10 100]./(fs/2), 'bandpass' );
lfp_filt = filtfilt(b_notch, a_notch, lfp);
lfp_filt = filtfilt(b_bp,a_bp, lfp_filt);

[b_notch,a_notch] = butter(4, [39 41]./(fs/2), 'stop' );
lfp_filt = filtfilt(b_notch, a_notch, lfp_filt);


[pow, f] = pmtm( lfp_filt, 2, nfft, fs);



subplot(2,1,1); hold on;
plot(lfp_filt, 'r', 'DisplayName', 'filtered'); 
legend('show', 'Location', 'eastoutside');

subplot(2,1,2); hold on;
plot(f, pow, 'r', 'DisplayName', 'filtered'); xlim([0 100]); 
legend('show', 'Location', 'eastoutside');


% clearvars pow f lfp b_notch a_notch b

%% 3. compare pairs of ex files regarding static and moving stimuli influence on the lfp

for i = 7:10%length(ex_all)/2
    expair = i;
    p = expair*2-1;
    
    figure('Name', [ex_all{p}.Header.fileName '-' ex_all{p+1}.Header.fileName]);
    sa(1) = subplot(2,2,1);
    lfp_dr_temp = lfpTimeDomain(ex_all{p}); title('drifting grating'); ylabel('\muV');
    lfp_dr(:, 1) = lfp_dr_temp(end, :);
    
    sb(1) = subplot(2,2,2);
    lfpFreqDomain(ex_all{p}, [1 100]); ylabel('Power \muV^2');
        
    sa(2) = subplot(2,2,3);
    lfp_st_temp = lfpTimeDomain(ex_all{p+1}); title('static grating'); xlabel('Time');
    lfp_st(:, 1) = lfp_st_temp(end, :);
    
    sb(2) = subplot(2,2,4);
    lfpFreqDomain(ex_all{p+1}, [1 100]); xlabel('Frequency');
    
    linkaxes(sa); linkaxes(sb);
    set(sb, 'YScale', 'log');

end


 



