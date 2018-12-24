function plot_encoding(serdata)
%%
% visualize the results from encoding analysis
% INPUT: serdata
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_drug.mat')
%

close all;
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_useufulfunc'))
addpath(genpath('Z:\Katsuhisa\code\integrated\cbrewer'))

% stimulus type
pairtype = serdata.session(459).results.pairtype;
row = [];
is5ht = [];
stmtype = [];
for i = 1:length(serdata.session)    
    if serdata.session(i).exist==1 && serdata.session(i).goodunit==1 ...
        && isstruct(serdata.session(i).results)
        stmtype = [stmtype, stm2val(serdata.session(i).results.stimulus)];
        row = [row, i];        
        if strcmp(serdata.session(i).results.drugname, '5HT')
            is5ht = [is5ht, 1];
        else
            is5ht = [is5ht, 0];
        end
    end
end

lenr = length(row);
disp(['Analized pairs of sessions: ' num2str(lenr)])
disp(['5HT sessions: ' num2str(sum(is5ht==1))])
disp(['NaCl sessions: ' num2str(sum(is5ht==0))])

serdata.session = serdata.session(row);

% pair type
switch pairtype
    case 'drug'
        names = {'baseline', 'drug'};
    case 'sc'
        names = {'low sc', 'high sc'};
    case 'ps'
        names = {'small ps', 'large ps'};
    case 'ps_drug'
        names = {'small ps drug', 'large ps drug'};
end

%%
% plot results from encoding analysis
stmcols = cbrewer('qual', 'Set2', 5);
switch pairtype
    case {'drug', 'ps_drug'}
        drugnames = {'NaCl', '5HT'};
        for s = 1:4
            for i = 1:2
                serdata_temp = serdata;
                serdata_temp.session = serdata_temp.session(is5ht==i-1 & stmtype==s);
                encode_visualizer(serdata_temp, names, stmcols(s,:))
                figname = [val2stm(s) '_' drugnames{i} '_' names{1} '_vs_' names{2} '_encoding'];
                set(gcf, 'Name', figname, 'NumberTitle', 'off')
                savefig(gcf, ['Z:\Katsuhisa\serotonin_project\LFP_project\' figname '.fig'])
            end
        end
    otherwise
        for s = 1:4
            serdata_temp = serdata;
            serdata_temp.session = serdata_temp.session(stmtype==s);
            encode_visualizer(serdata_temp, names, stmcols(s,:))
            figname = [val2stm(s) '_' names{1} '_vs_' names{2} '_encoding'];
            set(gcf, 'Name', figname, 'NumberTitle', 'off')
            savefig(gcf, ['Z:\Katsuhisa\serotonin_project\LFP_project\' figname '.fig'])
        end
end
        
% subfunction
function s = stm2val(stmtype)
switch stmtype
    case 'or'
        s = 1;
    case 'co'
        s = 2;
    case 'sf'
        s = 3;
    case 'sz'
        s = 4;
    case 'rc'
        s = 5;
end

function s = val2stm(val)
stms = {'or','co','sf','sz','rc'};
s = stms{val};

function encode_visualizer(serdata, names, col)
indnames = {'spike count', 'LFP (uV)', 'delta', 'theta', 'alpha', 'beta', 'gamma'};
fnames = {'effect','reliability', 'selectivity', 'snr2', 'discriminability', 'metabcost'};
leni = length(indnames);
lenf = length(fnames);
lenses = length(serdata.session);
% fz = 6;
figure;
c = 1;
for j = 1:lenf
    for i = 1:leni
        % scatter
        subplot(lenf, leni, c)
        v = nan(2, lenses);
        for k = 1:lenses            
            for l = 1:2
                try
                    if j == 1
                        v_temp = mean(serdata.session(k).results.cond(l).lfpstm.stm.tu{i}.mean);
                    else
                        v_temp = serdata.session(k).results.cond(l).lfpstm.stm.tu{i}.(fnames{j});
                        if length(v_temp) > 1
                            [~, prefidx] = max(serdata.session(k).results.cond(l).lfpstm.stm.tu{1}.mean);
                            v_temp = v_temp(prefidx);
                        end
                    end
                    if strcmp('selectivity', fnames{j})
                        v_temp = -log10(v_temp);
                    end
                catch
                    continue
                end
                v(l,k) = v_temp;
            end
        end
        nans = isnan(v(1,:)) | isnan(v(2,:));
        v(:, nans) = [];
        plot(v(1,:), v(2,:), 'o', 'color', col, 'markersize', 3)
        hold on;        
        % adjust axis
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
        xx = get(gca, 'XLim');
        yy = get(gca, 'YLim');
        minima = min([xx, yy]);
        maxima = max([xx, yy]);
        plot([minima maxima], [minima maxima], '-', 'color', 0.5*[1 1 1], ...
            'linewidth', 0.5)
        axis([minima maxima minima maxima])
        set(gca, 'XTick', [minima maxima], 'YTick', [minima maxima])
        if i==1
            if j==1
                xlabel(names{1})
                ylabel(names{2})
            else
                ylabel(fnames{j})
            end
        end        
        if j==1
            title(indnames{i})
        end
        % stats
        p = signrank(v(1,:), v(2,:));
        text(minima+0.1*(maxima-minima), maxima-0.1*(maxima-minima), ['p=' num2str(p)], ...
            'color', col)
        
        c = c + 1;
    end
end