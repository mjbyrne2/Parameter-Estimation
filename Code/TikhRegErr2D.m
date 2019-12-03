function e = TikhRegErr2D(dataSpec,operatorSpec,smoothingSpec,L,...
    origSpec)
% TikhRegErr2D determines the absolute error between origSpec, which is the 
% spectrum (DFT) of the unblurred test function, and the spectrum of
% the solution obtained through Tikhonov regularization with regularization
% parameter lambda (which is denoted solSpec below). 

% Initialization:
e = zeros(1,length(L));

for i = 1:length(e)
  
    regF_hat = conj(operatorSpec).*dataSpec./((abs(operatorSpec).^2) + ...
        ((L(i)^2)*smoothingSpec));
    regF = ifft2(regF_hat);
    F = ifft2(origSpec);

    % Relative error:
    e(i) = norm(regF-F)/norm(F);
    
end

end
