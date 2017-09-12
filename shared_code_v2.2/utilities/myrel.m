function e = myrel(V1,V2,p)

if nargin<3
    p = 2;
end

e = norm(col(V1-V2),p)/norm(col(V2),p);


