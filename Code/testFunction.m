function f = testFunction(dim,num)
% testFunction selects one of the test functions used in the numerical
% experiments.  The string dim specifies the dimension of the test 
% function ('1D', '2D', or '3D'), and the integer num specifies the test 
% function (1, 2, or 3). testFunction returns a row vector, not a function 
% handle.
%
% Companion files: periodicTestFunction_1D.m

if isequal(dim,'1D')
    
    N = 4096;
    t = linspace(0,1,N+1);
    t = t(1:end-1); % Equispaced N-point discretization of [0,1]
    
    if num == 1
        f = cos(4*pi*t).*sin(6*pi*t);
        return
    elseif num == 2
        f = periodicTestFunction_1D(N);
        return
    elseif num == 3
        f = cos(8*pi*t).*exp(sin(10*pi*t)-1);
        return
    end
    
elseif isequal(dim,'2D')
    
elseif isequal(dim,'3D')
    
end

end

