function [G_star] = make_star(n)

e1 = zeros(n,1);
e1(1) = 1;
G_star = e1 * ones(1,n) + ones(n,1) * e1.' - 2 * e1  * e1.';

end

