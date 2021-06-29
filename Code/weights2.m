function W = weights2(S,p,type)
% weights2.m generates weight matrices associated with windowed spectral
% Tikhonov regularization.
% - S is a matrix containing the (non-negative) spectra of an operator A
% - p is the an integer specifying the number of weight vectors
% - type specifies the type of weight vectors. Possibilites for type are
% 'linear', 'log', 'linearCosine', and 'logCosine'.

% Check number of input arguments:
switch nargin
    
    case 1  % Input is only S
        [m,n] = size(S);
        W = ones(m,n);
        return
        
    case 2  % No type specified
        if p == 1
            [m,n] = size(S);
            W = ones(m,n);
            return
        else
            type = 'linear';    % Default weight vectors are linear partitions
        end
        
    case 3
        if p == 1
            [m,n] = size(S);
            W = ones(m,n);
            return
        end
    
end

[m,n] = size(S);    % Dimension of S
W = zeros(m,n,p);   % 3D array of weight matrices
tolMax = 1/eps; % Tolerance used to measure the largest value in spectrum
if max(S(:)) > tolMax
    s = sort(S(:),'descend');   % Vectorize and sort S in descending order
    sMax = s(2);    % Second largest value in spectrum
else
    sMax = max(S(:));   % Largest value in spectum
end
sMin = max(min(S(:)),eps);   % Smallest value in spectrum (or eps)

% Check type of windows:
switch type
    case 'linear'
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(1) = sMax;    % Ensure the first partition value is sMax
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            W(:,:,k) = (omega(k) >= S) & (S > omega(k+1));
        end   
        
    case 'log'
        omega = logspace(log10(sMax),log10(sMin),p+1);  % Generate linear partitions
        omega(1) = sMax;    % Correct of numerical errors of using log10 with logspace
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            W(:,:,k) = (omega(k) >= S) & (S > omega(k+1));
        end     
        
    case 'linearCosine'
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(1) = sMax;    % Ensure the first partition value is sMax
        mid = (omega(1:end-1) + omega(2:end))/2;    % Determine midpoints
        omega(end) = -eps; % Set last partition value below zero
        
        % First and last windows:
        for j = 1:m
            for k = 1:n
                % First window:
                if (omega(1) >= S(j,k)) && (S(j,k) > mid(1))
                    W(j,k,1) = 1;
                elseif (mid(1) >= S(j,k)) && (S(j,k) > mid(2))
                    W(j,k,1) = cos((pi/2)*(mid(1)-S(j,k))/(mid(1)-mid(2)))^2;
                end
                % Last window:
                if (mid(end-1) >= S(j,k)) && (S(j,k) > mid(end))
                    W(j,k,p) = cos((pi/2)*(S(j,k)-mid(end))/(mid(end-1)-mid(end)))^2;
                elseif (mid(end) >= S(j,k)) && (S(j,k) > omega(end))
                    W(j,k,p) = 1;
                end
            end
        end
        
        if p > 2
            % Middle windows:
            for l = 2:p-1
                for j = 1:m
                    for k = 1:n
                        if (mid(l-1) >= S(j,k)) && (S(j,k) > mid(l))
                            W(j,k,l) = cos((pi/2)*(S(j,k)-mid(l))/(mid(l-1)-mid(l)))^2;
                        elseif (mid(l) >= S(j,k)) && (S(j,k) > mid(l+1))
                            W(j,k,l) = cos((pi/2)*(mid(l)-S(j,k))/(mid(l)-mid(l+1)))^2;
                        end  
                    end
                end
            end 
        end
        
    case 'logCosine'
        omega = logspace(log10(sMax),log10(sMin),p+1);  % Generate log partitions
        omega(1) = sMax;    % Correct of numerical errors of using log10 with logspace
%         mid = (omega(1:end-1) + omega(2:end))/2;    % Determine (linear) midpoints
        % Determine "log" midpoints:
        mid = zeros(1,p);
        for j = 1:p
            v = logspace(log10(omega(j)),log10(omega(j+1)),3);
            mid(j) = v(2);
        end
        omega(end) = -eps; % Set last partition value below zero
        
        % First and last windows:
        for j = 1:m
            for k = 1:n
                % First window:
                if (omega(1) >= S(j,k)) && (S(j,k) > mid(1))
                    W(j,k,1) = 1;
                elseif (mid(1) >= S(j,k)) && (S(j,k) > mid(2))
                    W(j,k,1) = cos((pi/2)*(mid(1)-S(j,k))/(mid(1)-mid(2)))^2;
                end
                % Last window:
                if (mid(end-1) >= S(j,k)) && (S(j,k) > mid(end))
                    W(j,k,p) = cos((pi/2)*(S(j,k)-mid(end))/(mid(end-1)-mid(end)))^2;
                elseif (mid(end) >= S(j,k)) && (S(j,k) > omega(end))
                    W(j,k,p) = 1;
                end
            end
        end
        
        if p > 2
            % Middle windows:
            for l = 2:p-1
                for j = 1:m
                    for k = 1:n
                        if (mid(l-1) >= S(j,k)) && (S(j,k) > mid(l))
                            W(j,k,l) = cos((pi/2)*(S(j,k)-mid(l))/(mid(l-1)-mid(l)))^2;
                        elseif (mid(l) >= S(j,k)) && (S(j,k) > mid(l+1))
                            W(j,k,l) = cos((pi/2)*(mid(l)-S(j,k))/(mid(l)-mid(l+1)))^2;
                        end  
                    end
                end
            end 
        end
                    
end

if max(S(:)) > tolMax
    W(1,1,1) = 1;    % Correct the weight of the largest component
end
      
end
