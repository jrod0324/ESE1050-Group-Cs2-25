function S = choose_leaders(G, P0, faulty, g, m)

%   Team: <Cs2-25>
%   Members: <Rohan, Jordan, Siddarth>
%   Strategy (2-3 sentences describing your rule, in words):
%   Want leaders that are the closest to the rally point, that way their
%   weights influence the follower robots to move more towards g. However,
%   we also want leaders to be ones with a lot of degrees, that way they
%   have more robots in their vicinity that can communicate and follow
%   them. As a result, we used both parameters via the sort function,
%   finding which points have the highest number of degrees, and then used
%   the smallest distance to the rally point to break ties. 

n = size(G, 1);
working = true(n, 1);
working(faulty) = false;

S = sortingCheck(G,P0,g,m,faulty);

% Calculates the number of degrees and distance to the rally point for each
% robot, return a list of m coordinate pairs that represent the leaders.
% Checks faulty indices by setting their degree values to -1, therefore
% they will never be called. 

function[orderedIndices] = sortingCheck(G,P0,g,m,faulty)
n = size(G, 1);
degrees = (G*ones(n,1));
degrees(faulty) = -1;
distances = vecnorm(P0 - g, 2,2);
robots = [distances,degrees]; 
[~,indices] = sortrows(robots,[2,1], {'descend','ascend'});
orderedIndices = indices(1:m);

end 
end
