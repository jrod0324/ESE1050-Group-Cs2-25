function G = make_disk(P,r)

    n = size(P,1);
    one = ones(n,1);

    s = sum(P.^2,2);

    delt = s*one' + one*s' - 2*(P*P');

    G = double(delt <= r^2);

    G(1:n+1:end) = 0;
 
end

Ptest = [0 0; 1 0; 2 0; 3 0; 5 0];
Gtest = make_disk(Ptest, 1.0)

