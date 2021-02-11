function x = xBigWin(alpha,W,D_hat,Delta,Lambda)

[~,~,p] = size(W);  % Number of windows
[ny,nx,r] = size(D_hat);  % Number of data sets
x = 0*D_hat;        % Initialization of x

switch p
    case 1  % p = 1  
        G = conj(Delta)./(repmat(abs(Delta).^2,1,1,r) + ...
            repmat((alpha*abs(Lambda).^2),1,1,r));
        x_hat = G.*D_hat;
        for l = 1:r
            x(:,:,l) = real(idct2(x_hat(:,:,l)));   % DCT of each matrix in the array
            % (dct2.m cannot handle 3D arrays yet, unlike fft2.m)
        end     
    otherwise  % p > 1
        I = sparse(repmat(eye(ny,nx),1,p)); % p-concatenation of identity matrices
        % Function that creates an (ny x nx x p) array of diagonal matrices using I:
        A = reshape(full(I*sparse(diag(reshape(repmat(alpha,ny,1),...
            ny*p,1)))),ny,nx,p);
        G = conj(Delta)./((abs(Delta).^2) + pagemtimes(abs(Lambda).^2,A));
        for l = 1:r
            x(:,:,l) = real(idct2(sum(G.*D_hat(:,:,l).*W,3)));
        end
end

end
