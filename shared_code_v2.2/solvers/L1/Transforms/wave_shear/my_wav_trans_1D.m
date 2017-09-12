function [D,Dt] = my_wav_trans_1D(wname,levels)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/22/2016

% these are only designed for 1-D problems
[Lo_d,Hi_d] = wfilters(wname);

D = @(x)wave_forward(x,Lo_d',Hi_d',levels);
Dt = @(c)wave_adjoint(c,Lo_d',Hi_d',levels);


function c = wave_forward(x,Lo_d,Hi_d,levels)

c = zeros(size(x,1),levels);
for i = 1:levels
    d = my_conv3D(x,Hi_d);
    x = my_conv3D(x,Lo_d);
    c(:,i) = d;
end



function x = wave_adjoint(c,Lo_d,Hi_d,levels)

c = reshape(c,numel(c)/levels,levels);
x = zeros(size(c));
for i = 1:levels
    x(:,i) = my_cconv3D(c(:,i),Hi_d);
    for j = i+1:levels
        c(:,j) = my_cconv3D(c(:,j),Lo_d);
    end
end
x = sum(x,2);



