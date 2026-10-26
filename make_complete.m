function [G_complete] = make_complete(n)

%   UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

G_complete = ones(n,1) * ones(1,n) - eye(n,n);

end
