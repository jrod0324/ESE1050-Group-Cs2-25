function G = make_disk(P, r)

s = sum(P.^2,2);
onus = ones(size(P, 1));
delta = s*onus' + onus*s' - 2*(P*P');

G = delta <= r.^2;

for i = 1:size(G, 1)
    for j = 1:size(G, 1)
        if i==j
            G(i, j) = 0;
        end
    end
end
end
