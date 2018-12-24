function exn = IRASA_single(ex, filename)
%%
% IRASA preprocessing for a single ex-file with preprocessed LFP data
%

% preset parameters
Fs = 1000; % sampling frequency
if ismember(1, contains(filename, 'xRC'))
    movingwin = [1000, 100]/Fs; % [window size, sliding step]
else
    movingwin = [250, 100]/Fs; % [window size, sliding step]
end
frange = [3 90];
win = movingwin(1)*Fs;
step = movingwin(2)*Fs;
folderName = '../../../LFP_project/Data/IRASAprepro/';
T = length(ex.Trials(end).LFP_prepro_time);
ntr = length(ex.Trials);

% apply IRASA
Frange = {[frange(1), 10], [15, 48]}; % define frequency range for power-law fitting
plawfnames = {'Beta', 'Cons', 'Plaw', 'Freq'};
nwin = floor((T - win)/step);
exn.IRASA_time = (0:step/Fs:step*nwin/Fs)';
for k = 1:ntr
    % separate fractal and oscillatory components using sliding window
    sig0 = zeros(win, nwin); sig1 = zeros(win, nwin);
    for i = 1:nwin+1
        sig0(:,i) = ex.Trials(k).LFP_prepro(...
            ceil((i-1)*step)+1:ceil((i-1)*step)+win);
        sig1(:,i) = ex.Trials(k).iLFP_prepro(...
            ceil((i-1)*step)+1:ceil((i-1)*step)+win);
    end
    Frac0_orig = amri_sig_fractal(sig0, Fs, 'detrend', 1, 'frange', frange);
    Frac1_orig = amri_sig_fractal(sig1, Fs, 'detrend', 1, 'frange', frange);
    
    % fitting power-law function to the fractal power spectra
    for j = 1:2
        Frac0 = amri_sig_plawfit(Frac0_orig, Frange{j});
        Frac1 = amri_sig_plawfit(Frac1_orig, Frange{j});
        for m = 1:4
            Frac0_orig.plawfit(j).(plawfnames{m}) = Frac0.(plawfnames{m});
            Frac1_orig.plawfit(j).(plawfnames{m}) = Frac1.(plawfnames{m});
        end
    end
    ex.Trials(k).IRASA = Frac0_orig;
    ex.Trials(k).iIRASA = Frac1_orig;
end
exn.Trials = ex.Trials;

% autosave =======================
filename = strrep(filename, '\', '/');
save([folderName filename], 'exn', '-v7.3')
disp([filename ' saved after IRASA!'])
