function x = xBigWin(alpha,p,D_hat,Delta,Lambda)

if p == 1
    x = BigGamma(alpha,p,D_hat,Delta,Lambda);
    for l = 1:size(x,3)
        x(:,:,l) = real(idct2(x(:,:,l)));   % DCT of each matrix in the array
        % dct2.m cannot handle 3D arrays yet, unlike fft2.m
    end
else
    
end

end