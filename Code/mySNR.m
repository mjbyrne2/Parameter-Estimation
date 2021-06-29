function S = mySNR(xReg,xTrue)
% mySNR calculates the SNR of a regularized solution xReg against the true
% solution xTrue. If xReg and xTrue are 3D arrays, mySNR calls the function
% arrayNorm.

r = size(xReg,3);   % Size of the third dimension

% Calculate SNR based on dimension:
switch r  
    case 1
        S = 20*log10(norm(xReg,'Fro')/norm(xReg-xTrue,'Fro'));      
    otherwise
        S = 20*log10(arrayNorm(xReg)./arrayNorm(xReg-xTrue));       
end

end
