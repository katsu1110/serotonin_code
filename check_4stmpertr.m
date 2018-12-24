function check_4stmpertr(ex)
% check pupil & lfp signals in the 4 stm/tr 
%

close all;
figure;
cols = jet(length(ex.Trials));
switch ex.exp.StimPerTrial
    case 1
        stmdur = 2;
    case 2
        stmdur = 0.45;
end
ntr = length(ex.Trials);
lenlfp = 10000;
lenps = 10000;
for i = 1:ntr
    if length(ex.Trials(i).LFP_prepro) < lenlfp
        lenlfp = length(ex.Trials(i).LFP_prepro);
    end
    if ex.Trials(i).Eye_prepro.enpos - ex.Trials(i).Eye_prepro.stpos + 1 < lenps
        lenps = ex.Trials(i).Eye_prepro.enpos - ex.Trials(i).Eye_prepro.stpos + 1;
    end
end
lfpv = nan(ntr, lenlfp);
psv = nan(ntr, lenps);
for i = 1:length(ex.Trials)
    % LFP
    subplot(2,1,1)
    t = ex.Trials(i).LFP_prepro_time + 0.5*(ex.Trials(i).labelseq - 1);
    y = ex.Trials(i).LFP_prepro;
    p = plot(t, y, '-', 'color', cols(i, :), 'linewidth', 1.5);
    p.Color(4) = 0.2;
    hold on;
    lfpv(i, :) = ex.Trials(i).LFP_prepro(end - lenlfp + 1:end);
    
    % pupil
    subplot(2,1,2)
    t = linspace(0, stmdur, ex.Trials(i).Eye_prepro.enpos - ex.Trials(i).Eye_prepro.stpos + 1); 
    t = t + 0.5*(ex.Trials(i).labelseq - 1);
    y = ex.Trials(i).Eye_prepro.ps(ex.Trials(i).Eye_prepro.stpos:ex.Trials(i).Eye_prepro.enpos);
    p = plot(t, y, '-', 'color', cols(i, :), 'linewidth', 3);
    p.Color(4) = 0.2;
    hold on;
    psv(i, :) = ex.Trials(i).Eye_prepro.ps(...
        ex.Trials(i).Eye_prepro.enpos - lenps + 1:ex.Trials(i).Eye_prepro.enpos);
end

% format
subplot(2,1,1)
ylabel('LFP')
xlim([0 2])
set(gca, 'box', 'off', 'tickdir', 'out')
subplot(2,1,2)
xlabel('time after stimulus onset (sec)')
xlim([0 2])
ylabel('pupil size')
set(gca, 'box', 'off', 'tickdir', 'out')
