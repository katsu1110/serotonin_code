function test_reconstruct_stLFP(data)
%%
% double-check whether stLFP can be reconstructed from 'data.mat'
% 
% data = convert_spiketimes(data);
close all;
figure;
for d = 1:2
    stlfp = getSTA(data{d}.calcium, 1:length(data{d}.calcium), data{d}.spike_times, 0.05, data{d}.fps);
    plot(-50:50, mean(stlfp, 1))
    hold on;
end

% function data = convert_spiketimes(data)
% for d = 1:2
%     data{d}.spikes = zeros(1, length(data{d}.calcium));
%     data{d}.spikes(data{d}.spike_times) = 1;
% end