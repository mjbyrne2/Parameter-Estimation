function G = BigGamma(alpha,p,D_hat,Delta,Lambda)
% Uses DCT

if p == 1
    G = conj(Delta).*D_hat./(repmat(abs(Delta).^2,1,1,size(D_hat,3)) + repmat((alpha*abs(Lambda).^2),1,1,size(D_hat,3)));
else
    
end

end
