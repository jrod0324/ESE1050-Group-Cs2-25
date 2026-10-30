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

% C.2 Spread and rendezvous time
function [delta, k_eps] = spread(X, Y)
n = size(X,1);
C = eye(n) - ones(n)/n;
delta = zeros(1, size(X,2));
for k = 1:size(X,2)
    delta(k) = sqrt(norm(C*X(:,k))^2 + norm(C*Y(:,k))^2);
end
k_eps = find(delta <= 0.01*delta(1), 1) - 1;
end









