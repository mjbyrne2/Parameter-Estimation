function x = myIFFT(x_hat)
% myIFFT.m is a rescaled version of the MATLAB built-in function ifft.m so
% that the transform corresponds with the matrix-vector multiplication of a
% vector with the unitary DFT matrix.

N = length(x_hat);
x = sqrt(N)*ifft(x_hat);

end

