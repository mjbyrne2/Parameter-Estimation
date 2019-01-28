function e = TikhRegErr(dataSpec,operatorSpec,smoothingSpec,L,...
    trunc,origSpec)
% TikhRegErr determines the absolute error between origSpec, which is the 
% spectrum (DFT) of the unblurred test function, and the spectrum of
% the solution obtained through Tikhonov regularization with regularization
% parameter lambda (which is denoted solSpec below). The nonnegative 
% integer trunc specifies the truncation level of solSpec.
%
% For this function to work properly, all spectrums should be zero-centered
% (see built-in fftshift.m).

% Initialization:
e = zeros(1,length(L));
N = length(dataSpec);
ind = (N-trunc)/2;  % Number of components to zero-out on both sides

for i = 1:length(e)
    
    filtCoeff = (conj(operatorSpec))./(abs(operatorSpec).^2 + ...
    (L(i).^2).*abs(smoothingSpec).^2);
    solSpec = filtCoeff.*dataSpec;

    % Truncate solution spectrum:
    solSpec([1:ind,N-ind:end]) = 0;

    % Absolute error:
    e(i) = err(solSpec,origSpec)*norm(origSpec);
    % e = sum((abs(original_spectrum - solution_spectrum)).^2); (original)

    % % Relative error:
    % e = err(solSpec,origSpec);
end

end
