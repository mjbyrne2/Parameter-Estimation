function D = dct1mtx(n)
%DCT2MTX Creates the DCT Type-I matrix

D = zeros(n);

for j = 1:n
    
    for k = 1:n
        
        if j == 1 || j == n
            a = 1/sqrt(2);
        else
            a = 1;
        end
        
        if k == 1 || k == n
            b = 1/sqrt(2);
        else
            b = 1;
        end
        
        D(j,k) = a*b*cos((pi/(n-1))*(j-1)*(k-1));
        
    end
end

D = sqrt(2/(n-1))*D;    % Scale the matrix

end
