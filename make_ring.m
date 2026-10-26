function [G_ring] = make_ring(n,k)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if k < 1 || k > n/2
    return
end

abs(eye(n) - ones(n))
G_ring = zeros(n);

e_j = zeros(n,1);
e_j(n) = 1;


S = zeros(n);

S(1,n) = 1; 
S(2:n,:) = eye(n-1,n);



for j = 1:k

    G_ring = G_ring + S + S';
    S = S * S;
end