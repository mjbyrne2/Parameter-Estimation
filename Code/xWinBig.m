function x = xWinBig(alpha,W,D_hat,delta,delta2,lambda)
% delta2 is modified from delta to avoid division by zero.

[~,~,p] = size(W);  % Number of windows
[ny,nx,r] = size(D_hat);  % Number of data sets
x = 0*D_hat;        % Initialization of x

switch p
    case 1  % p = 1  
        Phi = (delta.^2)./(repmat(delta.^2,1,1,r) + ...
            repmat((alpha*lambda.^2),1,1,r));
        x_hat = Phi.*D_hat./delta2;
        for l = 1:r
            x(:,:,l) = real(idct2(x_hat(:,:,l)));   % DCT of each matrix in the array
            % (dct2.m cannot handle 3D arrays yet, unlike fft2.m)
        end     
    otherwise  % p > 1
        I = sparse(repmat(eye(ny,nx),1,p)); % p-concatenation of identity matrices
        % Function that creates an (ny x nx x p) array of diagonal matrices using I:
        A = reshape(full(I*sparse(diag(reshape(repmat(alpha.^2,ny,1),...
            ny*p,1)))),ny,nx,p);    % A contains alpha^2
        Phi = (delta.^2)./((delta.^2) + pagemtimes(lambda.^2,A));
        x_hat = 0*x;    % Initalization of x_hat
        for l = 1:r
            x_hat(:,:,l) = sum(Phi.*D_hat(:,:,l).*W./delta2,3);
            x(:,:,l) = real(idct2(x_hat(:,:,l)));
        end
end

end
