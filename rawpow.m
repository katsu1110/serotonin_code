function [ x, y ] = rawpow( dat, varargin )

stimnorm = 0;
j=1;
while j<=length(varargin)
   switch varargin{j}
       case 'stimnorm'
           stimnorm = varargin{j+1};
   end
   j = j+1;
  
end



ind = [dat.is5HT];
% call the same function for baseline and 5ht/nacl data
pow_base = rawpow_helper(dat(ind), 'fname', stimnorm);
pow_drug = rawpow_helper(dat(ind), 'fname_drug', stimnorm);


f = 1:500/257:500;

if stimnorm
    for j = 1:size(pow_base, 3)
       
        subplot(size(pow_base, 3),2, j*2-1);
        plot(f, nanmean(pow_base(:,:,j), 1), 'r',...
            f, nanmean(pow_drug(:,:,j), 1), 'r--', 'LineWidth', 2);
        xlim([0, 100]);

    end
    
else
    plot(f, nanmean(pow_base, 1), 'r', f, nanmean(pow_drug, 1), 'r--', 'LineWidth', 2);
end

end




function pow = rawpow_helper(dat, filespec, stimnorm_flag)

if stimnorm_flag
    for unit = 1:length(dat)
        disp(['working on unit ' num2str(unit)]);
        fname = dat(unit).(filespec);
        ex = load(strrep(fname, dat(unit).cluster, 'lfp'), 'ex');
        ex = frequPreProc(ex.ex);
        for j= 1:length(dat(unit).ratepar)
            ind = dat(unit).ratepar(j) == [ex.Trials.(dat(unit).param1)];
            pow(unit, 1:257, j) = nanmean(horzcat(ex.Trials(ind).powstim), 2);
        end
    end
else
    parfor unit = 1:length(dat)
        disp(['working on unit ' num2str(unit)]);
        fname = dat(unit).(filespec);
        ex = load(strrep(fname, dat(unit).cluster, 'lfp'), 'ex');
        ex = frequPreProc(ex.ex);
        pow(unit, :, 1) = nanmean(horzcat(ex.Trials.powstim), 2);
    end
end


end


function ex = frequPreProc( ex )
% filters and plots data in the 

time = -0.5:0.001:0.5;

for ind = 1:length(ex.Trials)
    
    
    if isempty(ex.Trials(ind).LFP)
        ex.Trials(ind).powstim = nan(257, 1);
        ex.Trials(ind).powbase = nan(257, 1);
        
    else
        %% notch filter 
        Fs = 1000;          % sampling frequency
        wo = 50/(Fs/2);     % to be filtered frequency
        bw = wo/10;         % 50Hz / quality factor 

        [b,a] = iirnotch(wo,bw); 
        ex.Trials(ind).LFP_filt = filter(b, a, ex.Trials(ind).LFP);
        ex.Trials(ind).LFP_interp = ...
            interp1(ex.Trials(ind).LFP_ts, ex.Trials(ind).LFP_filt, time);

        
        %%  power spectrum
        x = ex.Trials(ind).LFP_interp(501:end);	%Sampled data in microvolt
        x2= ex.Trials(ind).LFP_interp(1:500);	%Sampled data in microvolt

        % mutli taper approach
        % power spectrum density returned by pmtm is normalized per unit frequency
        % in Chalk et al. time halfbandwidth was 3 and k was 5
        nw = 3;         % time half bandwidth product?????
        nfft = 2^nextpow2(length(x));

        ex.Trials(ind).powstim = pmtm(x, nw, nfft, Fs);
        ex.Trials(ind).powbase = pmtm(x2, nw, nfft, Fs);

    end
    
end


end







