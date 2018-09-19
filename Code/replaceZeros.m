function h = replaceZeros(x,r)
% replaceZeros returns the array h, where h is the array x except with all
% zeros in x replaced by r.

h = x;
h(h == 0) = r;

end
