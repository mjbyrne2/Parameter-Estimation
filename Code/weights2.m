function W = weights2(S,p,type)
% weights2.m generates weight matrices associated with windowed spectral
% Tikhonov regularization.
% - S is a matrix containing the (non-negative) spectra of an operator A
% - p is the an integer specifying the number of weight vectors
% - type specifies the type of weight vectors. Possibilites for type are
% 'linear', 'log', 'linearCosine', and 'logCosine'.

sMax = max(S(:));   % Largest value in spectum
sMin = max(min(S(:)),eps);   % Smallest value in spectrum (or eps)
[m,n] = size(S);    % Dimension of S
W = zeros(m,n,p);   % 3D array of weight matrices

% Check number of arguments:
if nargin == 2
    type = 'linear';    % Default weight vectors are linear partitions
end

% Check type of windows:
switch type
    case 'linear'
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            W(:,:,k) = (omega(k) >= S) & (S > omega(k+1));
        end   
        
    case 'log'
        omega = logspace(log10(sMax),log10(sMin),p+1);  % Generate linear partitions
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            W(:,:,k) = (omega(k) >= S) & (S > omega(k+1));
        end     
        
    case 'linearCosine'
        if p-2 < 0
            disp('p must be at least 2 to use cosine windows.')
            W = eye(m,n);
            return
        end
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(end) = -eps; % Set last partition value below zero
        mid = (omega(1:end-1) + omega(2:end))/2;    % Determine midpoints
        
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
        if p-2 < 0
            disp('p must be at least 2 to use cosine windows.')
            return
        end
        omega = logspace(log10(sMax),log10(sMin),p+1);  % Generate log partitions
        omega(end) = -eps; % Set last partition value below zero
        mid = (omega(1:end-1) + omega(2:end))/2;    % Determine midpoints
        
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
      
end
