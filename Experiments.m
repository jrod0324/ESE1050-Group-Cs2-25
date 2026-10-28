% experiments code

clear;
clc;
close all;

rng(105);
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];   
n = 100;                         
x0 = P0(:,1);         
y0 = P0(:,2);                   
r_min = 2.1;                    
K = 5000;

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



function W = make_WA(G)
n = size(G,1);
d = G*ones(n,1);                          
W = diag(1./(1+d)) * (eye(n) + G);        
end

function W = make_WB(G)
n = size(G,1);
d = G*ones(n,1);                          
L = diag(d) - G;                         
W = eye(n) - (1/max(d)) * L;              
end

function X = simulate_linear(W, x0, K)
X = zeros(length(x0), K+1);               
X(:,1) = x0;                              
for k = 1:K
    X(:,k+1) = W * X(:,k);               
end
if norm(X(:,end) - W^K*x0) > 1e-8 * norm(x0)   
    error('simulate_linear: mismatch');
end
end

function [delta, k_eps] = spread(X, Y)
n = size(X,1);
C = eye(n) - ones(n)/n;                   
delta = zeros(1, size(X,2));
for k = 1:size(X,2)
    delta(k) = sqrt(norm(C*X(:,k))^2 + norm(C*Y(:,k))^2);   
end
k_eps = find(delta <= 0.01*delta(1), 1) - 1;   
end



function [G_complete] = make_complete(n)
G_complete = ones(n,1) * ones(1,n) - eye(n,n);
end

function G = make_disk(P,r)
n = size(P,1);
one = ones(n,1);
s = sum(P.^2,2);
delt = s*one' + one*s' - 2*(P*P');
G = double(delt <= r^2);
G(1:n+1:end) = 0;
end



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


function [G_star] = make_star(n)

e1 = zeros(n,1);
e1(1) = 1;
G_star = e1 * ones(1,n) + ones(n,1) * e1.' - 2 * e1  * e1.';

end



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


