function plotGamma(ga1, x1, ga2, x2, fmt)

if nargin < 5
    fmt = '-ro';
end

m = length(x1);
n = length(x2);

hold on;

for i = 1:m
    plot([x1(i), x1(i)], [x2(1), x2(end)], '-b');
end

for j = 1:n
    plot([x1(1), x1(end)], [x2(j), x2(j)], '-b');
end

plot(ga1, ga2, fmt);

hold off;

axis([x1(1), x1(end), x2(1), x2(end)]);

end