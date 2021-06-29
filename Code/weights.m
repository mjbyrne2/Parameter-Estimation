function w = weights(s,p,type)
% weights.m generates weight vectors associated with windowed spectral
% Tikhonov regularization.
% - s is a vector containing the singular values of an operator A
% - p is the an integer specifying the number of weight vectors
% - type specifies the type of weight vectors. Possibilites for type are
% 'linear', 'log', 'linearCosine', and 'logCosine'.

n = length(s);  % Number of singular values
s  = sort(s,'descend'); % Sort singular values from largest to smallest
w = zeros(n,p); % Columns of w are the weight vectors to be returned
sMax = max(s);   % Largest singular value
sMin = max(min(s),eps);   % Smallest singular value (or eps)

% Check for Inf in s:
if isinf(sMax)
    disp('Error: Inf singular values detected.')
    return
end

% Check number of arguments:
if nargin == 2
    type = 'linear';    % Default weight vectors are linear partitions
end

% Check type of windows:
switch type
    case 'linear'
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(1) = omega(1) + eps;
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            w(:,k) = (omega(k) >= s) & (s > omega(k+1));
        end   
        
    case 'log'
        omega = logspace(log10(sMax),log10(sMin),p+1);  % Generate linear partitions
        omega(1) = omega(1) + eps;
        omega(end) = -eps; % Set last partition value below zero
        for k = 1:p
            w(:,k) = (omega(k) >= s) & (s > omega(k+1));
        end     
        
    case 'linearCosine'
        if p-2 < 0
            disp('p must be at least 2 to use cosine windows.')
            return
        end
        omega = linspace(sMax,sMin,p+1);  % Generate linear partitions
        omega(1) = omega(1) + eps;
        omega(end) = -eps; % Set last partition value below zero
        mid = (omega(1:end-1) + omega(2:end))/2;    % Determine midpoints
        
        % First and last windows:
        for j = 1:n
            % First window:
            if (omega(1) >= s(j)) && (s(j) > mid(1))
                w(j,1) = 1;
            elseif (mid(1) >= s(j)) && (s(j) > mid(2))
                w(j,1) = cos((pi/2)*(mid(1)-s(j))/(mid(1)-mid(2)))^2;
            end
            % Last window:
            if (mid(end-1) >= s(j)) && (s(j) > mid(end))
                w(j,p) = cos((pi/2)*(s(j)-mid(end))/(mid(end-1)-mid(end)))^2;
            elseif (mid(end) >= s(j)) && (s(j) > omega(end))
                 w(j,p) = 1;
            end
        end
        
        if p > 2
            % Middle windows:
            for k = 2:p-1
                for j = 1:n
                    if (mid(k-1) >= s(j)) && (s(j) > mid(k))
                        w(j,k) = cos((pi/2)*(s(j)-mid(k))/(mid(k-1)-mid(k)))^2;
                    elseif (mid(k) >= s(j)) && (s(j) > mid(k+1))
                        w(j,k) = cos((pi/2)*(mid(k)-s(j))/(mid(k)-mid(k+1)))^2;
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
        omega(1) = omega(1) + eps;
        omega(end) = -eps; % Set last partition value below zero
        mid = (omega(1:end-1) + omega(2:end))/2;    % Determine midpoints
        
        % First and last windows:
        for j = 1:n
            % First window:
            if (omega(1) >= s(j)) && (s(j) > mid(1))
                w(j,1) = 1;
            elseif (mid(1) >= s(j)) && (s(j) > mid(2))
                w(j,1) = cos((pi/2)*(mid(1)-s(j))/(mid(1)-mid(2)))^2;
            end
            % Last window:
            if (mid(end-1) >= s(j)) && (s(j) > mid(end))
                w(j,p) = cos((pi/2)*(s(j)-mid(end))/(mid(end-1)-mid(end)))^2;
            elseif (mid(end) >= s(j)) && (s(j) > omega(end))
                 w(j,p) = 1;
            end
        end
        
        if p > 2
            % Middle windows:
            for k = 2:p-1
                for j = 1:n
                    if (mid(k-1) >= s(j)) && (s(j) > mid(k))
                        w(j,k) = cos((pi/2)*(s(j)-mid(k))/(mid(k-1)-mid(k)))^2;
                    elseif (mid(k) >= s(j)) && (s(j) > mid(k+1))
                        w(j,k) = cos((pi/2)*(mid(k)-s(j))/(mid(k)-mid(k+1)))^2;
                    end             
                end
            end 
        end
                    
end
      
end
