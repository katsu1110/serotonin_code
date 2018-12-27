function psall(ex) 
close all
figure;

for i = 1:length(ex.Trials)
    subplot(1,2,1)
    plot(ex.Trials(i).Eye_prepro.psR(ex.Trials(i).Eye_prepro.stpos:ex.Trials(i).Eye_prepro.enpos))
    hold on;
    subplot(1,2,2)
    plot(ex.Trials(i).Eye_prepro.dpsR(ex.Trials(i).Eye_prepro.stpos:ex.Trials(i).Eye_prepro.enpos))
    hold on;
end