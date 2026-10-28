clear; clc; close all

rng(105);
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];
n = 100;


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


