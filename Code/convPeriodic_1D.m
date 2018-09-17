function xc = convPeriodic_1D(x,p)
% convPeriodic performs the one-dimensional periodic convolution of row 
% vector x and p (see details below). Vectors x and p are assumed to be row
% vectors of the same length N, with p being the zero-centered Gaussian 
% shifted from [-1/2,1/2] to [0,1] with the center now at 1/2. xc is the 
% resulting row vector of length N.
%
% The convolution is conducted in two equivalent ways. Using circular
% convolution, p is first used to construct h which is the the periodic 
% extension of the zero-centered Gaussian to the interval [0,1] (the 
% peaks are now at the ends of the vector).
% The second way is to use p directly but periodically extend x to the
% right and left of the original vector so that the periodic extension xp
% has length 2N. Next, the built-in conv.m is used with the 'valid' command
% to compute the linear convolution. Lastly, the last entry is trimmed so
% that the final result has length N. 
% As a note about the second method described previously, an equivalent way
% of using linear convolution is to use the command conv(p,xp,'same'),
% which automatically outputs a vector of length N. Michael, you should see
% if taking xp to be the 2N-1 extension of x (instead of 2N) and using the
% command conv(xp,p,'valid') would result in automatically outputing a 
% vector of length N.

% Circular convolution:
N = length(x);  % Should also be the length of p
h = [p((N/2)+1:end) p(1:N/2)];  % Periodic extension of zero-centered Gaussian
xc = cconv(x,h,N);

% % Linear convolution:
% xp = [x(end-floor(length(p)/2)+1:end) x x(1:floor(length(p)/2))];
% xc = conv(xp,p,'valid');
% xc = xc(1:length(x));
    
end