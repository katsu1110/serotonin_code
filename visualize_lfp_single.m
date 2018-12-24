function visualize_lfp_single(para)
%%
% plotting function to visualize the results from 'para_stmLFP.m' or
% 'para_stmLFP_mp.m'
%

close all;
fign = 1;

ss = length(para.stm.vals);
if ss==1
    cols = [0 0 0];
else
    cols = lines(ss);
end
if para.is5ht==1
    druglab = {'base', '5HT'};
else
    druglab = {'base', 'NaCl'};
end

% synchronization index
if isfield(para.cond(1), 'si')
    figure(fign);
    hcols = [0 0 0; 1 0 0];
    for s = 1:ss
        me = cell(1, 2);
        for d = 1:2
            subplot(1, ss, s)
            histogram(para.cond(d).si{s}, 'FaceColor', hcols(d, :));
            hold on;
            me{d} = nanmedian(para.cond(d).si{s});
        end
        yy = get(gca, 'YLim');
        for d = 1:2
            hold on;
            plot(me{d}*[1 1], yy, '-', 'color', hcols(d, :), 'linewidth', 2)
        end
        if s==1
            xlabel('trial-by-trial FF')
            ylabel('# trials')
        end
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
    set(gcf, 'Name', 'SI', 'NumberTitle', 'off')
    fign = fign + 1;
end

% LFP time-series
if isfield(para.cond(1), 'lfpfull')
    figure(fign);
    for s = 1:ss
        for d = 1:2
            % mean LFP
            subplot(1, 3, d)
            plot(para.ts, nanmean(para.cond(d).lfpfull{s}, 1), '-', 'color', cols(s, :))
            hold on;
        end
        % difference
        subplot(1, 3, 3)
        plot(para.ts, nanmean(para.cond(1).lfpfull{s}, 1) - nanmean(...
            para.cond(2).lfpfull{s}, 1), '-', 'color', cols(s, :))
        hold on;
    end 
    % format
    subplot(1,3,1)
    yy1 = get(gca, 'YLim');
    subplot(1,3,2)
    yy2 = get(gca, 'YLim');
    yy = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
    subplot(1,3,1)
    title(druglab{1})
    ylabel('LFP (uV)')
    xlim([-0.1 para.window{end}(2)])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    xx = get(gca, 'XLim');
    for s = 1:ss
        text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+(0.975 - 0.05*(s-1))*(yy(2)-yy(1)), ...
            num2str(para.stm.vals(s)), 'color', cols(s,:), 'fontsize', 6)
    end
    subplot(1,3,2)
    title(druglab{2})
    xlabel('time after stimulus onset (s)')
    xlim([-0.1 para.window{end}(2)])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(1,3,3)
    title('\Delta')
    yy = get(gca, 'YLim');
    xlim([-0.1 para.window{end}(2)])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    set(gcf, 'Name', 'LFP time-series', 'NumberTitle', 'off')
    fign = fign + 1;
end

% stLFP
if isfield(para.cond(1), 'sta')
    figure(fign);
    sta_t = linspace(-para.wnd, para.wnd, length(para.cond(1).sta.mean));
    spc_t = linspace(-para.wnd, para.wnd, length(para.cond(1).sta.t{end}));
    for s = 1:ss
        for d = 1:2
            % mean stLFP
            subplot(3, 3, d)
            plot(sta_t, para.cond(d).sta.mean(s, :), '-', 'color', cols(s, :))
            hold on;
            
%             % power spectrogram (freq x time)
%             if s==ss
%                 subplot(3, 3, d + 3)
%                 imagesc(spc_t, para.cond(d).sta.f{s}, para.cond(d).sta.p{s})
%             end
            
            % power (freq)
            subplot(3, 3, d + 6)
            plot(para.cond(d).sta.f{s}, nanmean(para.cond(d).sta.p{s}, 2), ...
                'color', cols(s, :))
            hold on;
        end
        % difference
        subplot(3, 3, 3)
        plot(sta_t, para.cond(1).sta.mean(s, :) - para.cond(2).sta.mean(s, :), '-', 'color', cols(s, :))
        hold on;
        subplot(3, 3, 6)
        try
            imagesc(spc_t, para.cond(1).sta.f{s}, para.cond(1).sta.p{ss} - para.cond(2).sta.p{ss})
        catch
        end
        subplot(3, 3, 9)
        plot(para.cond(d).sta.f{ss}, nanmean(para.cond(1).sta.p{ss}, 2) - ...
            nanmean(para.cond(2).sta.p{ss}, 2), 'color', cols(s, :))
        hold on;
    end 
    % format
    subplot(3,3,1)
    yy1 = get(gca, 'YLim');
    subplot(3,3,2)
    yy2 = get(gca, 'YLim');
    yy = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
    subplot(3,3,1)
    title(druglab{1})
    ylabel('LFP (uV)')
    xlim([-0.03 0.03])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    for s = 1:ss
        text(0.001, yy(1)+(0.975 - 0.08*(s-1))*(yy(2)-yy(1)), ...
            [num2str(para.stm.vals(s)) '(' num2str(para.cond(1).sta.nspk(s)) ')'], 'color', cols(s,:), 'fontsize', 6)
    end
    subplot(3,3,2)
    title(druglab{2})
    xlabel('time after spike (s)')
    xlim([-0.03 0.03])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    for s = 1:ss
        text(0.001, yy(1)+(0.975 - 0.08*(s-1))*(yy(2)-yy(1)), ...
            [num2str(para.stm.vals(s)) '(' num2str(para.cond(2).sta.nspk(s)) ')'], 'color', cols(s,:), 'fontsize', 6)
    end
    subplot(3,3,3)
    title('\Delta')
    yy = get(gca, 'YLim');
    xlim([-0.03 0.03])
    ylim(yy)
    hold on;
    plot([0 0], yy, '--', 'color', 0.5*[1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,4)
    ylabel('frequency (Hz)')
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,5)
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,6)
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,7)
    ylabel('power')
    set(gca, 'YScale', 'log')
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,8)
    xlabel('frequency (Hz)')
    set(gca, 'YScale', 'log')
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(3,3,9)
    set(gca, 'box', 'off', 'tickdir', 'out')
    colormap('jet')
    set(gcf, 'Name', 'stLFP', 'NumberTitle', 'off')
    fign = fign + 1;
end

% spectrogram
if isfield(para.cond(1), 'spectrogram')
    figure(fign);
    spc_t = linspace(-0.1, para.window{end}(2), length(para.cond(1).spectrogram.t{end}));
    for s = 1:ss
        for d = 1:2
            sz = size(para.cond(d).spectrogram.S{s});
            if sz(2) > sz(1)
                para.cond(d).spectrogram.S{s} = para.cond(d).spectrogram.S{s}';
            end
            
            % as it is... =========           
            
            % power (freq)
            subplot(4, 3, d)
            plot(para.cond(d).spectrogram.f{s}, 10*log10(nanmean(para.cond(d).spectrogram.S{s}(...
                para.window{end}(1) < spc_t & spc_t < para.window{end}(2), :), 1))', ...
                'color', cols(s, :))
            hold on;
            
            % power spectrogram (freq x time)
            if s==ss
                subplot(4, 3, d + 3)
                imagesc(spc_t, para.cond(d).spectrogram.f{s}, 10*log10(para.cond(d).spectrogram.S{s}'))
            end
            
            % response ==============
                       
            % power (freq)
            subplot(4, 3, d + 6)
            S = 10*log10(para.cond(d).spectrogram.S{s}');
            Sbase = nanmean(S(:, spc_t < 0), 2);
            Sres = S - repmat(Sbase, 1, length(spc_t));
            plot(para.cond(d).spectrogram.f{s}, nanmean(Sres(:, ...
                para.window{end}(1) < spc_t & spc_t < para.window{end}(2)), 2), ...
                'color', cols(s, :))
            hold on;
            
            % spectrogram
            if s==ss
                subplot(4, 3, d + 9)
                imagesc(spc_t, para.cond(d).spectrogram.f{s}, Sres)
            end
        end
        % difference
        subplot(4, 3, 3)
        plot(para.cond(d).spectrogram.f{s}, 10*log10(nanmean(para.cond(1).spectrogram.S{s}(...
            para.window{end}(1) < spc_t & spc_t < para.window{end}(2), :), 1)') - ...
            10*log10(nanmean(para.cond(2).spectrogram.S{s}(para.window{end}(1) < spc_t & ...
            spc_t < para.window{end}(2), :), 1)'), 'color', cols(s, :))
        hold on;
        subplot(4, 3, 6)
        if s == ss
            imagesc(spc_t, para.cond(d).spectrogram.f{s}, ...
                10*log10(para.cond(1).spectrogram.S{ss}') - 10*log10(para.cond(2).spectrogram.S{ss}'))
        end
        subplot(4, 3, 9)
        S = 10*log10(para.cond(1).spectrogram.S{s}');
        Sbase = nanmean(S(:, spc_t < 0), 2);
        Sres1 = S - repmat(Sbase, 1, length(spc_t));
        S = 10*log10(para.cond(2).spectrogram.S{s}');
        Sbase = nanmean(S(:, spc_t < 0), 2);
        Sres2 = S - repmat(Sbase, 1, length(spc_t));
        plot(para.cond(d).spectrogram.f{s}, nanmean(Sres1(:, ...
                para.window{end}(1) < spc_t & spc_t < para.window{end}(2)), 2) - ...
                nanmean(Sres2(:, para.window{end}(1) < spc_t & spc_t < para.window{end}(2)), 2), ...
                'color', cols(s, :))
        hold on;
        subplot(4, 3, 12)
        if s == ss
            imagesc(spc_t, para.cond(d).spectrogram.f{s}, Sres1 - Sres2)
        end
    end 
    % format
    for k = 1:12
        subplot(4, 3, k)
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
    subplot(4,3,1)
    yy1 = get(gca, 'YLim');
    subplot(4,3,2)
    yy2 = get(gca, 'YLim');
    yy = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
    subplot(4,3,1)
    title(druglab{1})
    ylim(yy)
    ylabel('power (dB)')
    xlabel('frequency (Hz)')
    hold on;
    subplot(4,3,2)
    title(druglab{2})
    ylim(yy)
    hold on;
    subplot(4,3,3)
    title('\Delta')
    subplot(4,3,4)
    xlabel('time after stimulus onset (sec)')
    ylabel('frequency (Hz)')
    title('overall power (dB)')
    subplot(4,3,7)
    ylabel('\Delta power (dB)')
    xlabel('frequency (Hz)')
    title('overall power (dB)')
    subplot(4,3,10)
    xlabel('time after stimulus onset (sec)')
    ylabel('frequency (Hz)')
    title('\Delta power (dB)')
    set(gcf, 'Name', 'spectrogram', 'NumberTitle', 'off')
    colormap('jet')
    fign = fign + 1;
end

% coherence
if isfield(para.cond(1), 'coherence')
    figure(fign);
    for s = 1:ss
        for d = 1:2
            % coherence
            subplot(2, 3, d)
            plot(para.cond(d).coherence.f{s}, para.cond(d).coherence.C{s}, '-', 'color', cols(s, :))
            hold on;
            
            % phase
            subplot(2, 3, d + 3)
            plot(para.cond(d).coherence.f{s}, para.cond(d).coherence.phi{s}, '-', 'color', cols(s, :))
            hold on;
        end
        % difference
        subplot(2, 3, 3)
        plot(para.cond(d).coherence.f{s}, para.cond(1).coherence.C{s} - ...
            para.cond(2).coherence.C{s}, '-', 'color', cols(s, :))
        hold on;
        subplot(2, 3, 6)
        plot(para.cond(d).coherence.f{s}, para.cond(1).coherence.phi{s} - ...
            para.cond(2).coherence.phi{s}, '-', 'color', cols(s, :))
        hold on;
    end 
    % format
    subplot(2,3,1)
    yy1 = get(gca, 'YLim');
    subplot(2,3,2)
    yy2 = get(gca, 'YLim');
    yyc = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
    subplot(2,3,4)
    yy1 = get(gca, 'YLim');
    subplot(2,3,5)
    yy2 = get(gca, 'YLim');
    yyp = [min([yy1(1), yy2(1)]), max([yy1(2), yy2(2)])];
    subplot(2,3,1)
    title(druglab{1})
    ylabel('coherence')
    ylim(yyc)
    hold on;
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2,3,2)
    title(druglab{2})
    ylim(yyc)
    hold on;
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2,3,4)
    ylim(yyp)
    ylabel('phase')
    xlabel('frequency (Hz)')
    hold on;
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2,3,5)
    ylim(yyp)
    hold on;
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2,3,3)
    title('\Delta')
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2,3,6)
    set(gca, 'box', 'off', 'tickdir', 'out')
    set(gcf, 'Name', 'coherence', 'NumberTitle', 'off')
    fign = fign + 1;
end