function x_hat = myFFT(x)

n = length(x);
x_hat = fft(x)/sqrt(n); % Rescale the fft to so it is now unitary

end

