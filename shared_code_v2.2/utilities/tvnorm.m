function tv = tvnorm(V,p)

if ~exist('p')
    p=1;
end

Dvx = [diff(V,1,2), V(:,1,:) - V(:,end,:)];
Dvy = [diff(V,1,1); V(1,:,:) - V(end,:,:)];
Dvz = cat(3, diff(V,1,3), V(:,:,1) - V(:,:,end));

tv = sum(sum(sum(abs(Dvx).^p))) + ...
        sum(sum(sum(abs(Dvy).^p))) + ...
        sum(sum(sum(abs(Dvz).^p)));
tv = tv^(1/p);
end