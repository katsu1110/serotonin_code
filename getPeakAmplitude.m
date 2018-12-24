function [a, t] = getPeakAmplitude( y, x )
% [a, t] = getPeak( x, y )
%
% compute the peak amplitude and the time of the peak occurance
%
% @CL 15.05.2017



% In case the signal can be flipped, the peak is a through.
% I did not consider this yet. 
[a, a_idx] = min(y); % find minimum == find negative peak 
t = x(a_idx); % time of peak


end

