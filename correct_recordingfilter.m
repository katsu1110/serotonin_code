function data = correct_recordingfilter(x, cutoff)
%%
% Correct filtering artifact by recording system
%
% NOTE that this function assumes that the recording system uses the 4th
% order butterworth filter (0-500Hz) to obtain LFP signals
%
% INPUT: x ... 1D vector LFP 
% OUTPUT: data ... corrected LFP
%
% https://github.com/open-ephys/analysis-tools/blob/master/lowFreqCorrection/loadAndCorrectPhase.m
%

% sampling rate (Hz)
Fs = 1000;

% % The values below come from direct measurement of the equipment, as described in the bioRxiv note
% switch cutoff
%   case 0.1
%     freq  = [0.0615    0.1214    0.1845    0.2247    0.2914    0.3732    0.8186    1.1102    2.0432    3.0051 11.1815   20.7900   30.1811];
%     phase = [3.0659    2.7502    2.4215    2.1768    2.0019    1.7454    1.1285    0.8774    0.5578    0.4007 0.1024    0.0389    0.0145];
%   case 1
%     freq  = [0.1067    0.1775    0.2308    0.2950    0.4092    0.8221    1.1241    2.0472    2.9354   11.2952  20.3804   29.6150];
%     phase = [3.6443    3.2698    2.8039    2.5902    2.2015    1.4175    1.1069    0.6644    0.4840   0.1121   0.0500    0.0213];
%   otherwise
%     error('loadAndCorrectPhase: this cutoff value isn''t supported')
% end
% freq(end+1) = 50;
% phase(end+1) = 0;

x = x - mean(x);
N = 2^nextpow2(min(numel(x), Fs*500)); % ~500s chunks

% Fourie Transform of the original signal
x_fft = fft(x);

% Fourie Transform of the butterworth filter
[b, a] = butter(4, 500/(Fs/2));


% freqX  = Fs * (1:N/2-1)/N; % frequency of each FFT component (to be computed below)
% freqXI = freqX >= freq(1) & freqX <= freq(end); % frequencies for which phase distortion was measured
% phaseDistort = zeros(N/2-1,1);
% phaseDistort(freqXI)  = interp1(log(freq), phase, log(freqX(freqXI)), 'pchip');
% phaseDistort(~freqXI) = 0; % The rationale is that for frequencies below freq(1) the power is so low it doesn't matter if we don't correct them.
% 
% data = zeros(size(x));
% for k = 1:ceil(numel(x)/N)
%   X = fft(x(N*(k-1)+1:min(N*k, numel(x))), N);    
%   X = X(2:N/2); 
%   X = abs(X).*exp(1i.*(angle(X) - phaseDistort)); % shift each frequency's phase by -phaseDistort
%   y = ifft([0; X; 0; conj(flipud(X))]);
%   data(N*(k-1)+1:min(N*k, numel(x))) = y(1:min(N*k, numel(x)) - N*(k-1));
% end