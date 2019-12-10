function y = myDCT(x)
% Computes the Type-II DCT of a vector x. The version of the Type-II DCT
% implemented is the one where the corresponding matrix is NOT orthogonal.
% This is in contrast to the MATLAB built-in function dct.m.

if isrow(x) == 1
    y = x;
    N = length(y);
    for j = 1:N
        y(j) = sum(x.*cos((pi/(2*N))*(((2*(1:N)))-1)*(j-1)));
    end
elseif iscolumn(x) == 1
    y = x;
    N = length(y);
    for j = 1:N
        y(j) = sum(x.*cos((pi/(2*N))*(((2*(1:N)'))-1)*(j-1)));
    end
else
    disp('Error: The input is not a row or column vector.')
end

end
