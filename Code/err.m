function e = err(D,T)
% err determines the relative error, calculated using the Euclidean norm,
% between the experimental data D and the true data T. T must be a row
% vector, though D can be a matrix with the same number of columns.
%
% If D and T are vectors, the error e is a scalar. If D is a matrix, then 
% the error is calculated for each row and the result e is a column
% vector. 

% Number of rows and columns of fE and fT:
[rowsD,colsD] = size(D);
[rowsT,colsT] = size(T);

% Test for the size of fE:
if isequal(size(D),size(T))
    e = norm(D - T)/norm(T);
    return
elseif colsD == rowsT
    T = T';
    Tmat = repmat(T,rowsD,1);  % Tmat is reshaped T
    e = sqrt(sum(abs(D-Tmat).^2,2));    % Euclidean norm of each row
    e = e./norm(T); % Relative error
    return
else
    Tmat = repmat(T,rowsD,1);  % Tmat is reshaped T
    e = sqrt(sum(abs(D-Tmat).^2,2));    % Euclidean norm of each row
    e = e./norm(T); % Relative error
    return
end
