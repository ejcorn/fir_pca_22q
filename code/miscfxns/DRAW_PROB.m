function F = DRAW_PROB(P,X,n)

C = cumsum(P);
F = nan(n,1);
for j = 1:n
    F(j) = X(1+sum(C(end)*rand>C));
end