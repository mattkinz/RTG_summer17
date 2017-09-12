function [overallerr, rel,relresids] = residualcomputation(angles,rec,stack)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


slices=size(stack,2);
if slices ~= size(rec,3)
    error('Reconstruction and stack appear different number of slices\n');
end

[W,~,bb]=getSIRT(angles,size(rec,1),stack);
residuals=zeros(slices,1);
nns = zeros(slices,1);
rec = reshape(rec,size(rec,1)*size(rec,2),slices);
err = W*rec-bb;
overallerr = norm(W*rec-bb,'fro')/norm(bb,'fro');

for i = 1:slices
    residuals(i)=norm(err(:,i));
    nns(i)=norm(bb(:,i));
end
relresids = residuals./(nns+1);
rel = median(relresids);
fprintf('Total relative error is %d\n',overallerr);

end
