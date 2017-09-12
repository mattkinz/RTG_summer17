function [V] = get_DT_weights(U,g)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

g = sort(g);
k = numel(g);
S = cell(k-1,1);
for i = 1:k-1
    S{i} = find(g(i)<=U & U<=g(i+1));
end

[a,b,c] = size(U);
U = U(:);
V = zeros(a*b*c,k);
V(S{1},1) = (U(S{1})-g(2))/(g(1)-g(2));
V(S{k-1},k) = (U(S{k-1})-g(k-1))/(g(k)-g(k-1));

for i = 2:k-1
    V(S{i-1},i) = (U(S{i-1})-g(i-1))/(g(i)-g(i-1));
    V(S{i},i) = (U(S{i}) - g(i+1))/(g(i) - g(i+1));
end

end