function Y = downsampleArray(X,rate)
% downsample2 is an extension of the built-in function downsample.m that
% downsamples an image (array) by the specified rate.

if isvector(X)
    Y = downsample(X,rate);
else
    Y = downsample(X,rate); % Downsample the columns
    Y = downsample(Y',rate);    % Downsample the rows
    Y = Y'; % Correct the orientation
end     

end
