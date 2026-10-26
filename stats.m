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
