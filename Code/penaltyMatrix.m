function Lambda = penaltyMatrix(name,n)
% penaltyMatrix.m returns spectrum Lambda of the penalty matrix specified 
% by the input name. Lambda is an n x n array computed using the 2D DCT
% which is used to blur an image through element-wise multiplication. If no
% value of n is given, the default value is 256. If no inputs are given,
% then Lambda is set to an array containing the names (as strings) of the
% possible penalty matrices.

validPenalties = ["Identity","Laplacian"];
nDefault = 256;

% Check inputs:
switch nargin
    case 0
        Lambda = validPenalties;
        return
    case 1
        if ~ismember(name,validPenalties)
            error('Error: The input penalty matrix is invalid.')
        else
            n = nDefault;   % Set default n
        end
    case 2
        if ~ismember(name,validPenalties)
            error('Error: The input penalty matrix is invalid.')
        end
    otherwise
        error('Error: Too many inputs in penaltyMatrix.m')
end

% Assign penalty matrix:
switch name
    case "Identity"
        Lambda = ones(n);  % DCT of l where L = I (Identity matrix)
    case "Laplacian" % Negative discrete Laplacian matrix
        L = zeros(n);
        cy = n/2;  % Row index of stencil center
        cx = n/2;  % Column index of stencil center
        L((cy-1):(cy+1),(cx-1):(cx+1)) = [0,-1,0;-1,4,-1;...
            0,-1,0];  % Place stencil within L
        e1 = zeros(n);
        e1(1,1) = 1;
        Lambda = dct2(dctshift(L,[cy,cx]))./dct2(e1);
end
        




end
