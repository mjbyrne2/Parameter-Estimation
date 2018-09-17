function f = periodicTestFunction_1D(N)
% periodicTestFunction_1D generates the second one-dimension test function,
% whose continuous version is periodic on the interval [0,1]. The input N
% specifies the number of points in the discretization and the output f is
% the discretized function in the form of a row vector. The function is
% defined piecewise below.

t = linspace(-pi,pi,N+1);
f = zeros(size(t));

for i = 1:N
    
    if -pi < t(i) && t(i) <= -pi/2
        f(i) = sin(4*(t(i)+pi));
    elseif -pi/2 < t(i) && t(i) <= -pi/3
        f(i) = 0;
    elseif -pi/3 < t(i) && t(i) <= -pi/4
        f(i) = 12*(t(i) + pi/3)/pi;
    elseif -pi/4 < t(i) && t(i) <= pi/4
        f(i) = 1;
    elseif pi/4 < t(i) && t(i) <= pi/3
        f(i) = -12*(t(i) - pi/3)/pi;
    elseif pi/3 < t(i) && t(i) <= pi/2
        f(i) = 0;
    elseif pi/2 < t(i) && t(i) <= pi
        f(i) = sin(4*(t(i)-pi/2));
    end

end

f = f(1:end-1);

end
