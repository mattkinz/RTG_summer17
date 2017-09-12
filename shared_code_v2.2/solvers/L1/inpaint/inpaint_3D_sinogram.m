function [U,out] = inpaint_3D_sinogram(bb,S,d,opts)



% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/22/2016



if numel(d)~=3
    error('this is 3D inpainting bro!');
end
%p = d(1);q = d(2); r = d(3);

if size(S,2)==3
    S = sub2ind(d,S(:,1),S(:,2),S(:,3));
end

if numel(S)~=numel(bb)
    error('number of data points and specified indices dont match');
end


A = @(x,mode)subdata_select(x,mode,S,d);
[opts.D,opts.Dt] = sparse_sino_trans(d,opts.W);
[U,out] = l1optimo(A,bb,d,opts);








function x = subdata_select(x,mode,S,d)

switch mode
    case 1
        x = x(S);
    case 2
        y = zeros(d);
        y(S) = x;
        x = y(:);
end


function [D,Dt] = sparse_sino_trans(d,W)

[F,Ft] = FD2D(3,d(1),d(2),1);
D = @(x)forward_t(x,W,F,d);
Dt = @(x)adjoint_t(x,W,Ft,d);


function y = forward_t(x,W,F,d)

%y = W'*x(:);
%y = x(:);
%y = col(F(y));
y = col(fft(reshape(x,d(1),d(2)),d(2),2));

function y = adjoint_t(x,W,Ft,d)

%y = Ft(x);
%y = W*col(y);
y = col(d(2)*ifft(reshape(x,d(1),d(2)),d(2),2));


