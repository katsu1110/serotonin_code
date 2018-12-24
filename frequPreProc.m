function ex = frequPreProc( ex )
% filters and plots data in the 

time = -0.5:0.001:0.5;



for ind = 1:length(ex.Trials)
    
    
    if isempty(ex.Trials(ind).LFP)
        ex.Trials(ind).powstim = nan(257, 1);
        ex.Trials(ind).powbase = nan(257, 1);
        ex.Trials(ind).powstim_f = 0;
        ex.Trials(ind).powbase_f = 0;
        
    else

        % notch filter 
        Fs = 1000;          % sampling frequency
        wo = 50/(Fs/2);     % to be filtered frequency
        bw = wo/10;         % 50Hz / quality factor 

        [b,a] = iirnotch(wo,bw); 
        ex.Trials(ind).LFP_filt = filter(b, a, ex.Trials(ind).LFP);
        ex.Trials(ind).LFP_interp = interp1(ex.Trials(ind).LFP_ts, ex.Trials(ind).LFP_filt, time);

        %%
        % fourier transformation and power spectrum

        x = ex.Trials(ind).LFP_interp(501:end);	%Sampled data in microvolt
        x2= ex.Trials(ind).LFP_interp(1:500);	%Sampled data in microvolt
        fs = 1000;
        n = length(x);
        fnq = fs/2 +1;

        nfft = 2^nextpow2(length(x));


        % normal fourier transformation
%         Y= fft(x, nfft);
%         Y = Y(1:nfft/2+1);
%         psdx = (1/(fs*length(x))) * abs(Y).^2;
%         psdx(2:end-1) = 2*psdx(2:end-1);
% 
%         f = 0:fs/nfft:fs/2;


        % mutli taper approach
        % power spectrum density returned by pmtm is normalized per unit frequency
        % in Chalk et al. time halfbandwidth was 3 and k was 5
        nw = 3;         % time half bandwidth product?????

        [ex.Trials(ind).powstim, ex.Trials(ind).powstim_f] = pmtm(x, nw, nfft, fs);
        [ex.Trials(ind).powbase, ex.Trials(ind).powbase_f] = pmtm(x2, nw, nfft, fs);

    end
    
end


end






