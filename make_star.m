function [G_ring] = make_star(n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
e1 = zeros(n,1);
e1(1) = 1;
G_ring = e1 * ones(1,n) + ones(n,1) * e1.' - 2 * e1  * e1.';

end
