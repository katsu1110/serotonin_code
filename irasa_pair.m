function [out0, out1] = irasa_pair(exns, stmtype, figoption)
%%
% IRASA analysis on the pair of exn-files
%

% input ====================
if nargin < 2; stmtype = 1; end
if nargin < 3; figoption = 0; end

% stimulus ==================
lasti = 3;
if stmtype == 1
    wnd = [0.8 2];
    stmtype = 'or';
else
    wnd = [0.2 0.45];
    stmtypes = {'or', 'or', 'co', 'sf', 'sz'};
    stmtype = stmtypes{stmtype};
end
stms = [exns{1}.Trials.(stmtype)];
unistm = unique(stms);
nuni = length(unistm);

% initialization =================
fields = {'IRASA', 'iIRASA'};
lenex = length(exns);
cols = jet(11);
lspc = {'-', '--'};
ts = exns{1}.IRASA_time;
ts = linspace(exns{1}.Trials(end).LFP_prepro_time(1), ...
    exns{1}.Trials(end).LFP_prepro_time(end), length(ts));
freq = exns{1}.Trials(end).(fields{1}).freq;
frac = cell(2, lenex);
osci = frac;
mixd = osci;
tumats = osci;
enc = tumats;
% bandnames = {'theta', 'alpha', 'beta', 'gamma'};
bandrange = {[3, 7], [8, 13], [14, 24], [25 48]};
lenb = length(bandrange);
pfitnames = {'Beta', 'Cons'};
paranames = {'fr', 'frac 3-7', 'frac 8-13', 'frac 14-24', 'frac 25-48', 'osci 3-7', 'osci 8-13', 'osci 14-24', 'osci 25-48'};
% lenp = length(paranames);

% data extraction ==================
lenf = length(freq);
lent = length(ts);
for f = 1:2
   for i = 1:lenex
      ntr = length(exns{i}.Trials);
      tumats{f, i} = nan(ntr, 2+2*lenb);
      frac{f, i} = zeros(lenf, lent, nuni);
      osci{f, i} = zeros(lenf, lent, nuni);
      mixd{f, i} = zeros(lenf, lent, nuni);
      stms = [exns{i}.Trials.(stmtype)];
      for n = 1:ntr
          % stimulus type
          idx = unistm == stms(n);
          
          % fractal & oscillation & mix
          if sum(isnan(exns{i}.Trials(n).(fields{f}).frac(:)))==0
                frac{f, i}(:, :, idx) = frac{f, i}(:, :, idx) + exns{i}.Trials(n).(fields{f}).frac;
          end
          if sum(isnan(exns{i}.Trials(n).(fields{f}).osci(:)))==0          
                osci{f, i}(:, :, idx) = osci{f, i}(:, :, idx) + exns{i}.Trials(n).(fields{f}).osci;      
          end
          if sum(isnan(exns{i}.Trials(n).(fields{f}).mixd(:)))==0                    
                mixd{f, i}(:, :, idx) = mixd{f, i}(:, :, idx) + exns{i}.Trials(n).(fields{f}).mixd;  
          end
          
          % tuning
          tumats{f, i}(n, 1) = stms(n); % stimulus type
          [~, tumats{f, i}(n, 2)]  = getSpks(exns{i}.Trials(n), [wnd(1), 0]); % firing rate
          
          % power law fit of fractal
          pos = 3;
          for k = 1:2
              for l = 1:2
                  tumats{f, i}(n, pos) = nanmean(...
                      exns{i}.Trials(n).(fields{f}).plawfit(k).(pfitnames{l})(end-lasti+1:end));
                  pos = pos + 1;
              end
          end
          
          % fractal & oscillation
%           pos = 3;
          f_tr = exns{i}.Trials(n).(fields{f}).freq;
          for b = 1:lenb
              frange = bandrange{b}(1) <= f_tr & f_tr <= bandrange{b}(2);
              tumats{f, i}(n, pos) = nanmean(nanmean(...
                   exns{i}.Trials(n).(fields{f}).frac(frange, end-lasti+1:end)));
              tumats{f, i}(n, pos+lenb) = nanmean(nanmean(...
                   exns{i}.Trials(n).(fields{f}).osci(frange, end-lasti+1:end)));
               pos = pos + 1;
          end          
      end
      
      % average fractal & oscillation & mix
      for s = 1:nuni          
          sntr = sum(stms == unistm(s));
          frac{f, i}(:, :, s) = frac{f, i}(:, :, s)/sntr;
          osci{f, i}(:, :, s) = osci{f, i}(:, :, s)/sntr;
          mixd{f, i}(:, :, s) = mixd{f, i}(:, :, s)/sntr;
      end
      
      % tuning analysis
      if wnd(2) - wnd(1) < 1
          tumat = tumats{f, i};
          tumat(tumat(:,1) >= 1000, :) = [];
          for b = 1:size(tumat, 2) - 1
              try
                enc{f, i}{b} = encoding_tuning(...
                        tumat(:, 1), tumat(:, 1+b), stmtype);
              catch
                  continue
              end
          end
      end
   end
   
  % visualize ===============================
  if figoption == 1     
      figure(f);
      if wnd(2) - wnd(1) > 1          
          % visualize fractal & oscillation
          for i = 1:2
              for s = 1:nuni
                  subplot(1, 2, 1)
                  m = squeeze(mixd{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', 'b')
                  hold on; 
                  m = squeeze(frac{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', 'r')
                  hold on;
                  subplot(1, 2, 2)
                  m = squeeze(mixd{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', 'b')
                  hold on; 
                  m = squeeze(frac{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', 'r')
                  hold on;
                  m = squeeze(osci{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', [0 1 0]/2)
                  hold on;
              end
          end    
          
          % format
         subplot(1, 2, 1)           
         xlabel('frequency (Hz)')
         ylabel('mixed')
         set(gca, 'XScale', 'log', 'YScale', 'log')         
         set(gca, 'box', 'off', 'tickdir', 'out')
         xx = get(gca, 'XLim');
         yy = get(gca, 'YLim');
         text(xx(1)+0.45*(xx(2)-xx(1)), yy(1)+0.95*(yy(2)-yy(1)), 'total', 'color', 'b')
         text(xx(1)+0.45*(xx(2)-xx(1)), yy(1)+0.6*(yy(2)-yy(1)), 'fractal', 'color', 'r')
         subplot(1, 2, 2)          
         xlabel('frequency (Hz)')
         ylabel('oscillation')          
         set(gca, 'box', 'off', 'tickdir', 'out')
         xx = get(gca, 'XLim');
         yy = get(gca, 'YLim');
         text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.95*(yy(2)-yy(1)), 'total', 'color', 'b')
         text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.85*(yy(2)-yy(1)), 'fractal', 'color', 'r')
         text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+0.75*(yy(2)-yy(1)), 'oscillation', 'color', [0 1 0]/2)
      else          
          if strcmp(stmtype, 'co')
              sidx = zeros(1, 3);
              [~, sidx(1)] = min(abs(unistm - 0.25));
              [~, sidx(2)] = min(abs(unistm - 0.5));
              [~, sidx(3)] = min(abs(unistm - 1));
              lcols = cols([1 6 11], :);
          else
              lcols = cols([1 11], :);
              mfr = zeros(1, nuni);
              for s = 1:nuni          
                  mfr(s) = nanmean(tumats{f, i}(tumats{f, i}(:,1)==unistm(s), 2));
              end
              [~, sortidx] = sort(mfr);
                sidx = [sortidx(1), sortidx(end)];
                if sortidx(1)==1
                    sidx(1) = sortidx(2);
                end
          end
          lens = length(sidx);
          
          % visualize fractal & oscillation
          lenp = size(tumat, 2) - 1;
          for i = 1:2
              for s = 1:lens
                  subplot(2, lenp, 1:6)
                  m = squeeze(mixd{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', lcols(s,:))
                  hold on;
                  m = squeeze(frac{f, i}(:, :, s));
                  plot(freq, mean(m(:, end-lasti+1:end), 2), lspc{i}, 'color', lcols(s,:))
                  hold on;
                  subplot(2, lenp, 7:lenp)
                  m = squeeze(osci{f, i}(:, :, s));
                  plot(freq, mean(m(:, ts >= wnd(1) & ts <= wnd(2)), 2), lspc{i}, 'color', lcols(s,:))
                  hold on;
              end
          end

          % format
          subplot(2, lenp, 1:6)          
         xlabel('frequency (Hz)')
         ylabel('mixed')
         set(gca, 'XScale', 'log', 'YScale', 'log')         
         set(gca, 'box', 'off', 'tickdir', 'out')
         subplot(2, lenp, 7:lenp)          
         xlabel('frequency (Hz)')
         ylabel('oscillation')          
         set(gca, 'box', 'off', 'tickdir', 'out')
         
         % visualize tuning curve
         tu_cols = [0 0 0; 1 0 0];
        for d = 1:lenp
            for i = 1:2
                subplot(2, lenp, d + lenp)
                errorbar(enc{f, i}{d}.unistm, enc{f, i}{d}.mean, ...
                    enc{f, i}{d}.std./sqrt(enc{f, i}{d}.ntr), 'color', tu_cols(i, :), 'capsize', 0);
                hold on;
            end   
            title(paranames{d})
         end
      end        
      set(gcf, 'Name', fields{f}, 'NumberTitle', 'off')
  end
end

% outputs ==============================
out0.ts = ts;
out0.freq = freq;
out0.frac = {frac{1, 1}, frac{1, 2}};
out0.osci = {osci{1, 1}, osci{1, 2}};
out0.mixd = {mixd{1, 1}, mixd{1, 2}};
out0.tumats = {tumats{1,1}, tumats{1, 2}};
out0.enc = {enc{1,1}, enc{1, 2}};
out1.ts = ts;
out1.freq = freq;
out1.frac = {frac{2,1}, frac{2, 2}};
out1.osci = {osci{2,1}, osci{2, 2}};
out1.mixd = {mixd{2,1}, mixd{2, 2}};
out1.tumats = {tumats{2,1}, tumats{2, 2}};
out1.enc = {enc{2,1}, enc{2, 2}};