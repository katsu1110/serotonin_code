function  LFPplotBatch( exSpkin, exLFPin, exSpkin2, exLFPin2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% settings
nw = 3;     % number of tapers in the multitaper approach
params.Fs = 1000; params.trialave = 1;
params.tapers = [nw nw*2-1];
params.fpass = [1 100];
params.time = -0.05:1/params.Fs:getstimdur(exSpkin);

% filter settings
notchfreq = [48 52];    notchorder = 2; % notch filter for 5Hz frequency
lowpfilter = 250;       lowporder = 2;  % general low pass filter


% remove empty rows
[exLFPin, exSpkin] =  rmemptyrows(exLFPin, exSpkin);
[exLFPin2, exSpkin2] =  rmemptyrows(exLFPin2, exSpkin2);

exlfp = frequAnalysis(exLFPin, 'notchf', notchfreq, 'notchord', notchorder, ...
    'lowpf', lowpfilter, 'lowpord', lowporder, 'nw', nw);
exlfp2 = frequAnalysis(exLFPin2, 'notchf', notchfreq, 'notchord', notchorder, ...
    'lowpf', lowpfilter, 'lowpord', lowporder, 'nw', nw);

plotLFP(exlfp, exlfp2, params)
plotLFP_spkavg(exSpkin, exlfp, exSpkin2, exlfp2)
end

%%

function plotLFP(exlfp, exlfp_drug, params)
% LFP in the time and frequency domain
figure('tag', 'lfp', 'Position', [100 296   962   435])
subplot(4,2,1); lfpave_base = lfpTimeDomain(exlfp); title('notch and lowpass filtered LFP'); legend('off')
subplot(4,2,3); lfpave_drug = lfpTimeDomain(exlfp_drug); legend('off')
subplot(4,2,[5 7]); plot(params.time, nanmean(lfpave_base), getCol(exlfp_drug)); hold on;
    plot(params.time, nanmean(lfpave_drug), getCol(exlfp_drug), 'LineStyle', '--');
    xlabel('time (s)'); xlim([params.time(1) params.time(end)]); 
legend('baseline', 'drug', 'Location', 'NorthOutside', 'Orientation', 'horizontal'); 
 crossl;
 
 
fidx = exlfp.Trials(1).FREQ>= params.fpass(1) & exlfp.Trials(1).FREQ <= params.fpass(2);
f = exlfp.Trials(1).FREQ(fidx);

subplot(4,2,2); powave_base = lfpFreqDomain(exlfp, params.fpass); title('Power plot')
subplot(4,2,4); powave_drug = lfpFreqDomain(exlfp_drug, params.fpass);  legend('off')
subplot(4,2,[6 8]); plot(f, nanmean(powave_base), getCol(exlfp_drug)); hold on;
    plot(f, nanmean(powave_drug), getCol(exlfp_drug), 'LineStyle', '--');
set(gca,'YScale', 'log'); xlabel('Frequency (Hz)'); xlim(params.fpass);
legend('baseline', 'drug', 'Location', 'NorthOutside', 'Orientation', 'horizontal'); 
 crossl;
end

%%
function plotLFP_spkavg(exspk, exlfp, exspk_drug, exlfp_drug)
% spike averaged LFP
figure('tag', 'spk_avg', 'Position', [581   314   689   659])


s1 = subplot(4,1,1); spktriglfp( exspk, exlfp, 'plot', true, 'rawflag', true);
spktriglfp( exspk, exlfp, 'plot', true); ylim auto
legend('show', 'Location', 'EastOutside');     crossl;


s2 = subplot(4,1,2); spktriglfp( exspk_drug, exlfp_drug, 'plot', true, 'rawflag', true);
spktriglfp( exspk_drug, exlfp_drug, 'plot', true); ylim auto
set(findobj(s2, 'type', 'line'), 'Color', getCol(exspk_drug));
legend('show', 'Location', 'EastOutside');     crossl;

                
s3 = subplot(4,1,3); 
spktriglfp( exspk, exlfp, 'plot', true);hold on
spktriglfp( exspk_drug, exlfp_drug, 'plot', true);
s3.Children(1).Color = getCol(exspk_drug);
ylim auto ;    crossl;


s4 = subplot(4,1,4); 
spktriglfp( exspk, exlfp, 'plot', true, 'rawflag', true);hold on
spktriglfp( exspk_drug, exlfp_drug, 'plot', true, 'rawflag', true);
s4.Children(1).Color = getCol(exspk_drug); ylim auto;     crossl;


end

function plotLFP_coh(exspk, exlfp, exspk_drug, exlfp_drug, params)
% spike LFP coherence plot

subplot(3,1,1); [coh_base, fcoh_base] = spkfieldcoh(exspk, exlfp, params); xlim(params.fpass);
subplot(3,1,2); [coh_drug, fcoh_drug] = spkfieldcoh(exspk_drug, exlfp_drug, params); xlim(params.fpass);

subplot(3,1,3);
spkfieldcoh(exspk, exlfp, params); hold on;
set(findobj(ax_spkfieldcoh_cmp,'type', 'line'), 'Color', 'k');
spkfieldcoh(exspk_drug, exlfp_drug, params);
xlim(params.fpass);
legend('baseline', 'drug', 'Location', 'EastOutside');
end
%%

function [exLFPin, exSpkin] =  rmemptyrows(exLFPin, exSpkin)
% remove all trials with empty lfp files
idx = ~ cellfun(@isempty, {exLFPin.Trials.LFP});
exLFPin.Trials = exLFPin.Trials(idx); exSpkin.Trials = exSpkin.Trials(idx);
end

function [stim, vals] = getStimVals(exSpkin,exSpkin2)
% stimulus parameters and corresponding values
[stim, valsB] = getStimParam(exSpkin);
[~, valsD] = getStimParam(exSpkin2);
vals = intersect(valsD, valsB); % find the common stimuli 
end



function stimdur = getstimdur(ex)

fname = getFname(ex);

if isempty(strfind(fname, 'RC'))
    stimdur = 0.45;
else
    stimdur = 2;
end

end

