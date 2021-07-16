function v = arrayNorm(A)
% For a 3D array A, arrayNorm returns a row vector v consisting of the
% Frobenius norms of A(:,:,r) for r = 1:size(A,3).

R = size(A,3);  % Size of third dimension
switch R
    case 1  % If A is a matrix and not a 3D array
        v = norm(A,'fro');
    otherwise
        v = zeros(1,R);
        for r = 1:R
            v(r) = norm(A(:,:,r),'fro');
        end
end

end

