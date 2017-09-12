function [b,scl] = Scaleb(b)

% Scales mu and f so that the finite difference of f is neither too small 
% nor too large.
%
% If option is assigned, mu will be scaled accordingly.
%
% Written by: Chengbo Li
[m,n,k] = size(b);
b = b(:);
threshold1 = .5;      % threshold is chosen by experience.
threshold2 = 1.5;
scl = 1;
b_dif = abs(max(b) - min(b));

if b_dif < threshold1
    scl = threshold1/b_dif;
    b = scl*b;
else if b_dif > threshold2
        scl = threshold2/b_dif;
        b = scl*b;
    end
end
b = reshape(b,m,n,k);
return