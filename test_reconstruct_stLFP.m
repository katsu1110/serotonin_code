function test_reconstruct_stLFP(data)
%%
% double-check whether stLFP can be reconstructed from 'data.mat'
% 
% data = convert_spiketimes(data);
close all;
figure;
col = {'k', 'r'};
stlfp = cell(1,2);
for i = 1:length(data)
    stlfp_tr = getSTA(data{i}.calcium, 1:length(data{i}.calcium), data{i}.spike_times, 0.05, data{i}.fps);
    if data{i}.cell_num == 1
        stlfp{1} = [stlfp{1}; stlfp_tr]; 
    else
        stlfp{2} = [stlfp{2}; stlfp_tr];
    end
end
for c = 1:2
    plot(-50:50, mean(stlfp{c}, 1), 'color', col{c})
    hold on;
end

% function data = convert_spiketimes(data)
% for d = 1:2
%     data{d}.spikes = zeros(1, length(data{d}.calcium));
%     data{d}.spikes(data{d}.spike_times) = 1;
% end