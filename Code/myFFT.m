function x_hat = myFFT(x)
% myFFT.m is a rescaled version of the MATLAB built-in function fft.m so
% that the transform corresponds with the matrix-vector multiplication of a
% vector with the unitary DFT matrix.

N = length(x);
x_hat = fft(x)/sqrt(N);

end

