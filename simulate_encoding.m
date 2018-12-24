function simulate_encoding(varargin)
%%
% simulation encoding property using membrane potential model
%

addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc\neuro'))

% generate spikes using membrane potential model
v_thre = [1 3 6 8 11];
n_thre = length(v_thre);
disp('generating spike activity')
mm = mem2fire('stmtype', 'or', 'ntr', 200, 'v_thre', v_thre);

% format matrix (stm, res) & encoding analysis
disp('encoding analysis...')
tu = cell(1, n_thre);
for n = 1:n_thre
    % matrix formation
    begin = mm.stm.ntr*(n-1) + 1;
    mat = [mean(mm.stm.stm, 2), ...
       mean(mm.res.spk(begin:begin + mm.stm.ntr -1, :), 2)];
   
    % encoding analysis
    tu{n} = encoding_tuning(mat(:, 1), mat(:, 2));
end

% visualize
close all;
figure;
cols = lines(n_thre);
indnames = {'reliability', 'selectivity', 'snr2', 'discriminability', 'metabcost'};
lenin = length(indnames);
for n = 1:n_thre
    % tuning curve
    subplot(2,3,1)
    hold on;
    errorbar(tu{n}.unistm, tu{n}.mean, tu{n}.std, ...
        '-', 'color', cols(n,:), 'capsize', 0)
    if n==n_thre
       ylabel('spike response')
       xlabel('stimulus')
       set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
    end
    
    % indices
    for i = 1:lenin
        subplot(2,3,1+i)
        hold on;
        v = tu{n}.(indnames{i});
        if length(v) > 1
           prefidx = tu{n}.mean==max(tu{n}.mean);
           v = v(prefidx);
        end
        if strcmp(indnames{i}, 'selectivity')
            v = -log(v);
        end
        hold on;
        bar(n, v, 'facecolor', cols(n,:), 'edgecolor', 'w')
        if n==n_thre
           ylabel(indnames{i})
           if i==lenin-1
              xlabel('Vm threshold')              
           end
           set(gca, 'XTick', 1:n_thre, 'XTickLabel', v_thre)
           set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); 
        end
    end
end

