%% All code for Parts 1, 2
clear; clc; close all

%% Tests on small 5-robot examples
% Expected: complete E=10 diam=1 | star E=4 diam=2 | ring(k=1) E=5 diam=2
test_complete = graph_stats(make_complete(5))
test_star     = graph_stats(make_star(5))
test_ring     = graph_stats(make_ring(5,1))
% Robots on a line at x = 0,1,2,3,5 with r = 1 -> E=3, diam=Inf
test_disk = graph_stats(make_disk([0 0; 1 0; 2 0; 3 0; 5 0], 1))

% Protocol test on 5-robot star: A should meet at x = 20/13, B at x = 2
Gt = make_star(5);
xt = [0; 1; 2; 3; 4];
XA_test = simulate_linear(make_WA(Gt), xt, 50);
XB_test = simulate_linear(make_WB(Gt), xt, 50);
final_A = XA_test(:,end)'
final_B = XB_test(:,end)'

%% Part 1: Communication Networks
%% B.2 Deployment
rng(105);
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];
n  = 100;          % total number of robots
x0 = P0(:,1);      % starting x
y0 = P0(:,2);      % starting y
K  = 5000;         % simulation steps

%% B.3 r_min search
r = 0.1;
stats = graph_stats(make_disk(P0, r));
while stats.diam == inf
    r = r + 0.1;
    stats = graph_stats(make_disk(P0, r));
end
r_min = round(r, 1)

%% B.3 Build the four networks
G1 = make_complete(n);
G2 = make_star(n);
G4 = make_disk(P0, r_min);
stats_G4 = graph_stats(G4);
k  = round(stats_G4.E / n)
G3 = make_ring(n, k);
stats_G1 = graph_stats(G1);
stats_G2 = graph_stats(G2);
stats_G3 = graph_stats(G3);

%% B.4 Table deliverable
fprintf('\n%-10s %8s %10s %10s %10s\n', 'Network', 'Edges', 'MeanDeg', 'MaxDeg', 'Diameter')
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Complete', stats_G1.E, mean(stats_G1.d), max(stats_G1.d), stats_G1.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Star',     stats_G2.E, mean(stats_G2.d), max(stats_G2.d), stats_G2.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Ring',     stats_G3.E, mean(stats_G3.d), max(stats_G3.d), stats_G3.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Disk',     stats_G4.E, mean(stats_G4.d), max(stats_G4.d), stats_G4.diam)

%% B.4 T4 network plot
figure
plot(graph(G4), 'XData', P0(:,1), 'YData', P0(:,2))
axis equal
xlabel('x'); ylabel('y')
title(['T4 disk network, r = ' num2str(r_min)])

%% Part 2: Rendezvous Protocols
%% C.2 Experiments
G5 = make_disk(P0, r_min + 0.4);
Gs = {G1, G2, G3, G4, G5};
names = {'Complete','Star','Ring','Disk r_{min}','Disk r_{min}+0.4'};

% Experiment 1: How fast?
figure; hold on
for t = 1:5
    W = make_WA(Gs{t});
    X = simulate_linear(W, x0, K);
    Y = simulate_linear(W, y0, K);
    [delta, k_eps] = spread(X, Y);
    semilogy(0:K, delta)
    fprintf('%-18s k_eps = %d\n', names{t}, k_eps)
end
set(gca, 'YScale', 'log')
xlabel('step k'); ylabel('\delta[k]')
legend(names); title('E1: Protocol A spread')

% Experiment 2: Where do they meet?
for t = 1:4
    G = Gs{t};
    d = G*ones(n,1);
    v = 1 + d;
    predA = [v'*x0, v'*y0] / sum(v);   % weighted average
    predB = [mean(x0), mean(y0)];      % plain average
    WA = make_WA(G);
    WB = make_WB(G);
    XA = simulate_linear(WA, x0, K);  YA = simulate_linear(WA, y0, K);
    XB = simulate_linear(WB, x0, K);  YB = simulate_linear(WB, y0, K);
    [~, kA] = spread(XA, YA);
    [~, kB] = spread(XB, YB);
    fprintf('\n%s\n', names{t})
    fprintf('  A: predicted (%.3f, %.3f), simulated (%.3f, %.3f), k_eps = %d\n', ...
        predA, XA(1,end), YA(1,end), kA)
    fprintf('  B: predicted (%.3f, %.3f), simulated (%.3f, %.3f), k_eps = %d\n', ...
        predB, XB(1,end), YB(1,end), kB)
end

% Star: how far from the average does Protocol A meet?
d = G2*ones(n,1);
v = 1 + d;
predA_star = [v'*x0, v'*y0] / sum(v);
dist_from_avg = norm(predA_star - [mean(x0), mean(y0)])

% D1: Steady State for simulate_anchored
D  = make_deployment(1);
n  = size(D.G,1);  d = D.G*ones(n,1);
W  = diag(1./(1+d))*(eye(n)+D.G);
x0 = D.P0(:,1);  y0 = D.P0(:,2);  g = D.g(:);

lam = ones(n,1);  lam(D.faulty) = 0;
xs = steady_state(W, lam, x0);
ys = steady_state(W, lam, y0);
X  = simulate_anchored(W, lam, x0, x0, 5000);
Y  = simulate_anchored(W, lam, y0, y0, 5000);
err = max(norm(xs - X(:,end)), norm(ys - Y(:,end)));
fprintf(mean(xs), mean(ys), g(1), g(2))

%% ================= Functions (must be at the end) =================
% T1 Broadcast (complete)
function G = make_complete(n)
G = ones(n,1) * ones(1,n) - eye(n);
end

% T2 Base station (star)
function G = make_star(n)
e1 = zeros(n,1);
e1(1) = 1;
G = e1 * ones(1,n) + ones(n,1) * e1' - 2 * e1 * e1';
end

% T3 Ring
function G = make_ring(n,k)
if k < 1 || 2*k >= n
    error('make_ring: need 1 <= k < n/2');
end
S = zeros(n);
S(1,n) = 1;
S(2:n,:) = eye(n-1,n);     % S*e_j = e_(j+1)
G = zeros(n);
Sj = eye(n);
for j = 1:k
    Sj = Sj * S;           % S^j
    G = G + Sj + Sj';
end
end

% T4 Radio range (disk graph)
function G = make_disk(P,r)
n = size(P,1);
one = ones(n,1);
s = sum(P.^2,2);
delt = s*one' + one*s' - 2*(P*P');
G = double(delt <= r^2);
G(1:n+1:end) = 0;
end

% B.4 Structure from matrix products
function stats = graph_stats(G)
n = size(G,1);
one = ones(n,1);
I = eye(n);
stats.d = G*one;
stats.E = (1/2)*one'*G*one;
A = I + G;
R = A;
stats.diam = inf;
for k = 1:n
    if all(R(:) > 0)
        stats.diam = k;
        return
    end
    R = double(A*R > 0);
end
end

% C.1 Protocol A
function W = make_WA(G)
n = size(G,1);
d = G*ones(n,1);
W = diag(1./(1+d)) * (eye(n) + G);
end

% C.1 Protocol B
function W = make_WB(G)
n = size(G,1);
d = G*ones(n,1);
L = diag(d) - G;
W = eye(n) - (1/max(d)) * L;
end

% C.1 Simulator

% C1 Experiment estimates: 

% Estimate T1: Since we are removing the whole sets average using the
% broadcast technique should drop off quickly to zero and then become
% constant

% Estimate T2: Again all of the robots are able to communicate, but
% communication has to go through more channels, which means it will
% converge to zero quickly but not as quickly as T1

% Estimate T3: First case where all robots are not able to hear each other,
% only their k nearest neighbors. As all the robots become closer, they
% will gain more and more neighbors and the value will again converge to
% zero semi quickly, slower than the first two cases.

% Estimate T4 r min: This is the best case scenario, should del(k) should
% seem to decrease linearly when shown on a log scale graph

% Estimate T4 r min +0.4 : This should show similar behavior to r min, with
% a slightly different intercept
% 

function X = simulate_linear(W, x0, K)
X = zeros(length(x0), K+1);
X(:,1) = x0;
for k = 1:K
    X(:,k+1) = W * X(:,k);
end
if norm(X(:,end) - W^K*x0) > 1e-8 * norm(x0)
    error('simulate_linear: last column does not match W^K*x0');
end
end


% C1 Experiment Results: 

% T1: The prediction for T1 is correct. All robots are able to see the
% average and should move to the average in one step with del(k) being zero
% with some error, which is exactly what happened. The gap between the two
% zones has little effect

% T2: The star method also converged to the average very quickly but
% slightly slower than T1 as predicted. In the space that the star
% algorithm builds, each robot should be no more than 2 communications away
% from the next. For every step, the stars will become move halfway towards
% the center and become more connected. With this specific outline, a
% specific robots number of links does not matter, it is the spread of info
% across the hub that determines convergence speed. Again the gap between
% zones has little effect

% T3: The logic for T3 was slightly incorrect. Each robot is only allowed
% to 2k neighbors, so information is spread slower with each time step
% which means more steps are needed to converge. Additionally, each of
% these steps is a local average so it takes more time for robots on
% opposite sides of the setup to communicate with each other on movements.
% Again, the gap between zones has little effect on the final result

% T4 r min: This is the slowest convergence rate even though it has the
% same amount of edges as a ring. This is mainly because only a handful of
% rings cross the gap, which creates a bottleneck in information. I was
% correct that it would be slightly slower than using a larger radius value

%T4 r min +0.4 : I was slightly wrong on this front because the larger
%radius creates more links for information to travel across different
%zones, so it does not exhibit the same slow convergence as the other T4
%case

% C.2 Spread and rendezvous time

% Rendezvous point for each protocol: PA creates a weighted sum, so the
% swarms meets at the calculated degree weighted average of starting
% positions, with weight being based on d_i for each robot i. In contrast,
% protocol uses the ones vector which creates an ordinary average of
% initial positions. In only part A, robots with more neighbors will have
% more influence on where the swarm converges to.

% Predicted vs simulated meeting points are shown in command windows output

% Star average meeting point: Protocol will meet at the average exactly,
% because it is just a plain sum

% Protocol A has the hub(lead robot) has the most weight(1+99) while all
% other robots have weight(1+1 = 2). With 100 robots the hub carries around
% 33% of total weight and each leaf only carries 0.66%. This means the
% meeting point is pulled 33% of the way towards the initial robot. So you
% can use the initial robots place and the average vector to find the
% actual conversion point for Protocol A

% Protocol reaches a meeting point in 7 steps while protocol B takes 450
% steps. A is much faster because each leaf only has one neighbor(the
% initial robot) and moves to a midpoint between that and its initial
% position. This is able to average all positions with 1 step. Meanwhile
% with the total average, the total distance covered by one robot towards
% the hub is 1/99, for the 99 non leading robots. This smaller distance
% covered per step helps to create the difference in convergence speeds
% between the two protocols


function [delta, k_eps] = spread(X, Y)
n = size(X,1);
C = eye(n) - ones(n)/n;
delta = zeros(1, size(X,2));
for k = 1:size(X,2)
    delta(k) = sqrt(norm(C*X(:,k))^2 + norm(C*Y(:,k))^2);
end
k_eps = find(delta <= 0.01*delta(1), 1) - 1;
end


% D1 Leaders and Faulty Robots
function X2 = simulate_anchored(W, lambda, a, x0, K)

greekA = diag(lambda);
X2 = zeros(length(x0), K+1);
X2(:, 1) = x0;

for i = 1:K
    X2(:,i+1) = (greekA*W*X2(:,i)) + (eye(length(x0)) - greekA)*a;
end

end

function ss = steady_state(W, lambda, a);
Lam = diag(lambda);
ss  = (eye(length(a)) - Lam*W) \ ((eye(length(a)) - Lam)*a);
end


%D2 The challenge
function S = choose_leaders(G, P0, faulty, g, m)

%   Strategy (2-3 sentences describing your rule, in words):
%   We want leaders that are the closest to the rally point, that way their
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

% D3: Scoring
function J = score_mission(D, S)
m = 4;                  % num leaders
T = 125;                % num steps
n = size(D.G, 1);       % number of robots
g = D.g;                % rally point

%leader set check
if length(S) > m
    error('Too many leaders');
end
for i = 1:length(S)
    for j = 1:length(D.faulty)
        if S(i) == D.faulty(j)
            error('A leader is faulty');
        end
    end
    for j = i+1:length(S)
        if S(i) == S(j)
            error('Leader repitition');
        end
    end
end


working = ones(n,1);
working(D.faulty) = 0;

% Setup
n = size(D.G,1);
d = (D.G)*ones(n,1);
W = diag(1./(1+d)) * (eye(n) + D.G);
x0 = D.P0(:,1);
y0 = D.P0(:,2);

% No leaders test, establish e_base for cost calculation
lam = ones(n,1);        % followers
lam(D.faulty) = 0;      % faulty robots don't move
a = x0;                 
b = y0;
X = simulate_anchored(W, lam, a, x0, T);
Y = simulate_anchored(W, lam, b, y0, T);
e_base = zeros(1,T);
for k = 1:T
    dist2 = (X(:,k+1) - g(1)).^2 + (Y(:,k+1) - g(2)).^2;  % euclidian distance of robots
    e_base(k) = sum(working .* dist2) / sum(working);      
end

% With leaders
lam(S) = 0.5;           % leaders
a(S) = g(1);            % leaders anchored at g
b(S) = g(2);
X = simulate_anchored(W, lam, a, x0, T);
Y = simulate_anchored(W, lam, b, y0, T);
e_new = zeros(1,T);
for k = 1:T
    dist2 = (X(:,k+1) - g(1)).^2 + (Y(:,k+1) - g(2)).^2;
    e_new(k) = sum(working .* dist2) / sum(working);
end

% Score
J = (100/T) * sum(1 - e_new ./ e_base);
end
