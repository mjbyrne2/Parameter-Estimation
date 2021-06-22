function xMin = fmingrid(f,x0,pts,iter)
% fmingrid find an approximation to a local minimizer of the scalar-valued
% function f within a domain whose edges are specified by x0. The dimension
% of x0 is n x 2, where n is the dimension of the input of f. f is a
% function of a single argument but this argument can be vector of length
% n. pts is the maximum number of points used to approximate the domain and
% iter is the number of iterations of the minimization process.

[n,~] = size(x0);  % Dimension of input of f
nx = floor(nthroot(pts,n));  % Size of discretization of each dimension
N = nx^n;   % Total number of domain points
F = zeros(N,1); % Intialized function evaluation
C = cell(1,n);  % Initialize cell array used to create discretized domain

% Loop over specified number of iterations:
for i = 1:iter
    
    % Construct the cell containing the individual domain discretizations:
    for d = 1:n
    %     D = logspace(log10(x0(d,1)),log10(x0(d,2)),nx);   % Log discretization
        D = linspace(x0(d,1),x0(d,2),nx);  % Linear discetization
        D(1) = x0(d,1);    % Correct first value
        D(end) = x0(d,2);  % Correct last value
        C{d} = D;
    end

    % Convert from cell to double:
    [G{1:n}] = ndgrid(C{:});    % Construct a grid-like array
    for k = n:-1:1
        X(:,k) = G{k}(:);
    end
    
    % Evaluate the function at each point:
    for j = 1:N
        F(j) = f(X(j,:));   
    end

    % Find the minimizer:
    indMin = find(F == min(F)); % Row index of the minimizer
    xMin = X(indMin,:);

    % Redefine x0 for use in any further iterations:
    for d = 1:n       
        % Check if minimum occurs on the boundary:
        if sum(x0(d,:) == xMin(d)) > 0
            disp('Minimum detected on the boundary. Minimization has ceased.')
            return
        end
        L = min(abs(x0 - repmat(xMin(:),1,2)),[],'all');    % Find the minimum distance to the boundary
        x0 = repmat(xMin(:),1,2) + repmat([-L,L],n,1);  % Redefine x0    
    end
    
    % Clear G and X for any further iterations:
    clear G X
    
end

end
