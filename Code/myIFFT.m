function x = myIFFT(x_hat)

n = length(x_hat);
x = sqrt(n)*ifft(x_hat); % Rescale the ifft to so it is now unitary

end

