function y = ir_dct2(x)
% todo: comments etc

[nx, ny, np] = size(x);
y = permute(reshape(dct(reshape(permute(reshape(dct(reshape(x,nx,[])),nx,ny,np),[2,1,3]),ny,[])),ny,nx,np),[2 1 3]);

