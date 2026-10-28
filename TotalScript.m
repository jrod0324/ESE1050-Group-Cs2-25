%% All code for Parts 1,2
clear; clc; close all

% Part 1: Communication Networks

%% B.2 Deployment: 
rng(105);
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];
n = 100;  % total number of robots



%% B.3 Four Networks

% T1 Broadcast(complete)
function [G_complete] = make_complete(n)

G_complete = ones(n,1) * ones(1,n) - eye(n,n);

end

% T2 Base station(star)

function [G_star] = make_star(n)

e1 = zeros(n,1);
e1(1) = 1;
G_star = e1 * ones(1,n) + ones(n,1) * e1.' - 2 * e1 * e1.';

end

% T3 Ring

function [G_ring] = make_ring(n,k)
    if k < 1 || k > n/2
    return
    end

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
end

% T4 Radio Range(Disk Graph)
function G = make_disk(P,r)
    n = size(P,1);
    one = ones(n,1);
    s = sum(P.^2,2);
    delt = s*one' + one*s' - 2*(P*P');
    G = double(delt <= r^2);
    G(1:n+1:end) = 0;
end

% r_min optimization
r = 0.1;
stats = graph_stats(make_disk(P0, r));
while stats.diam == inf
    r = r + 0.1;
    stats = graph_stats(make_disk(P0, r));
end
r_min = round(r, 1);

G1 = make_complete(n);
G2 = make_star(n);
G4 = make_disk(P0, r_min);
stats_G4 = graph_stats(G4);
k = round(stats_G4.E / n) ;    
G3 = make_ring(n, k);

stats_G1 = graph_stats(G1);
stats_G2 = graph_stats(G2);
stats_G3 = graph_stats(G3);

% table deliverable
fprintf('\n%-10s %8s %10s %10s %10s\n', 'Network', 'Edges', 'MeanDeg', 'MaxDeg', 'Diameter')
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Complete', stats_G1.E, mean(stats_G1.d), max(stats_G1.d), stats_G1.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Star',     stats_G2.E, mean(stats_G2.d), max(stats_G2.d), stats_G2.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Ring',     stats_G3.E, mean(stats_G3.d), max(stats_G3.d), stats_G3.diam)
fprintf('%-10s %8d %10.2f %10d %10d\n', 'Disk',     stats_G4.E, mean(stats_G4.d), max(stats_G4.d), stats_G4.diam)

% T4 network plotted
figure
plot(graph(G4), 'XData', P0(:,1), 'YData', P0(:,2))
axis equal
xlabel('x'); ylabel('y')
title(['T4 disk network, r = ' num2str(r_min)])


%% B.4 Structure from Matrix Products
function stats = graph_stats(G)

n = size(G,1);
one = ones(n,1);
I = eye(n);

deg_vec = G*one;
edge_count = (1/2)*(one')*(G)*(one);

stats.d = deg_vec;
stats.E = edge_count;

A = I + G;
R = A;
stats.diam = inf;

for k = 1:n
    if all(R(:) > 0)    
        stats.diam = k;
        return
    end
    R = double(A * R > 0);   
end
end



%% C.1 Two protocols

function X = simulate_linear(W, x0, K)


X = zeros(length(x0), K+1); % X dimensions
X(:,1) = x0; % First col of X is x0
% setting the next columns of X through the matrix vector product
for k = 1:K
    X(:,k+1) = W * X(:,k);
end

% checks last column of X is correct
if norm(X(: , end) - W^K * x0) > 0
    error('Error');
end
end

D = diag(d);
L = D - G;
I = eye(n);

WA = diag(1./(1+D)) * (I+G);
e = 1/(max(d));
WB = I - e*L;

XA = simulate_linear(WA, x0, K);
XB = simulate_linear(WB, x0, K);

%% C.2 Experiments

G1 = make_complete(n);
G2 = make_star(n);
G4 = make_disk(P0, r_min);
k = round(graph_stats(G4).E/n);
G3 = make_ring(n,k);
G5 = make_disk(P0, r_min +0.4);

Gs = {G1, G2, G3, G4, G5};
names = {'Complete','Star','Ring','Disk r_{min}','Disk r_{min}+0.4'};


% Experiment 1

figure;
hold on

for t = 1:5
    W = make_WA(Gs{t});
    X = simulate_linear(W, x0, K);
    Y = simulate_linear(W ,y0, K);
    [delta, k_eps] = spread(X, Y);
    semilogy(0:K, delta) % plot spread(ensure its log scale)
end

set(gca, 'YScale', 'log')
xlabel('step k');
ylabel('delta[k]');
legend(names);
title('E1: Protocol A spread')

% Experiment 2
for t = 1:4
    G = Gs{t};
    d = G*ones(n,1);

    v = 1+d;
    predA = [v'*x0, v'*y0]/sum(v);
    predB = [mean(x0), mean(y0)];

    WA = make_WA(G);
    WB = make_WB(G);
    XA = simulate_linear(WA, x0, K);
    YA = simulate_linear(WA, y0, K);
    XB = simulate_linear(WB, x0, K);
    YB = simulate_linear(WB, y0, K);
    [~, kA] = spread(XA,YA);
    [~, kB] = spread(XB,YB);

    fprintf('\n%s\n', names{t})
    fprintf('  A: predicted (%.3f, %.3f), simulated (%.3f, %.3f), k_eps = %d\n', ...
        predA, XA(1,end), YA(1,end), kA)
    fprintf('  B: predicted (%.3f, %.3f), simulated (%.3f, %.3f), k_eps = %d\n', ...
        predB, XB(1,end), YB(1,end), kB)
end




