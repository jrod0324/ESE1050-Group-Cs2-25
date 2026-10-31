function S = choose_leaders(G, P0, faulty, g, m)
%CHOOSE_LEADERS  Choose which robots become leaders (ESE 105, Case Study II, Part III).
%
%   S = choose_leaders(G, P0, faulty, g, m)
%
%   Inputs
%     G       n-by-n adjacency matrix of the communication network
%     P0      n-by-2 landing positions; row i is robot i's (x, y)
%     faulty  vector of indices of robots with broken drives
%     g       1-by-2 rally point
%     m       maximum number of leaders
%
%   Output
%     S       vector of at most m distinct robot indices, none of them faulty
%
%   Rules (see Section D.4 of the handout)
%     - Use only the network's structure: anything you can compute from G,
%       P0, faulty and g with the tools from Part I (degrees, hop counts,
%       distances, products and powers of G, ...).
%     - Do NOT simulate the swarm, and do NOT call score_mission.
%     - Must run in under 2 seconds on one deployment.
%     - Do not change the function name or the order of the inputs.
%       Helper functions are allowed, but put them in this file, below
%       this function, so that this single file is all we need to run it.
%
%   Team: <team name>
%   Members: <Rohan, Jordan, Siddarth>
%
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

% ---- replace this placeholder with your strategy ----


k = false;
while (k == false)
ranking = sortingCheck(G,P0,g,m);
    if all(ranking) == 1
        k = true;
    end
end

S = ranking;   
end



function[orderedIndices] = sortingCheck(G,P0,g,m)
n = size(G, 1);
degrees = (G*ones(n,1));
distances = vecnorm(P0 - g, 2,2);

robots = [distances,degrees]; 

[~,indices] = sortrows(robots,[2,1], {'descend','ascend'});

orderedIndices = indices(1:m);

end 
