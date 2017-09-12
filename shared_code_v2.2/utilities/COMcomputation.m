function c = COMcomputation(volume)

[m,n,k]= size(volume);

M = sum(sum(sum(volume)));

s1 = sum(sum(volume,3),2);

s1 = (1:m)*s1/M;

s2 = sum(sum(volume,3),1);

s2 = sum(s2.*(1:n)/M);


s3 = sum(sum(volume,2),1);

s3 = sum(reshape(s3,1,k).*(1:k)/M);


c = [s1 s2 s3];

end