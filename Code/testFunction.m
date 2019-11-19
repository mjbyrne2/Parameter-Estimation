function f = testFunction(num)
% testFunction(num) selects one of the 1D test functions used in the 
% numerical experiments.  The integer num specifies the 
% test function (1, 2, or 3). testFunction returns a function handle.
%
% Companion files: periodicTestFunction_1D.m

if num == 1
    f = @(t) cos(4*pi*t).*sin(6*pi*t);
    return
elseif num == 2
    f = @(t) periodicTestFunction_1D(N);
    return
elseif num == 3
    f = @(t) cos(8*pi*t).*exp(sin(10*pi*t)-1);
    return
end

end

