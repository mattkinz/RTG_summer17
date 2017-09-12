function [rec] = multigrid_rec(stack,angles,opt)

opt = check_multi_opt(opt);

[bb,W,v1,v2,slices] = get_multigrid(stack,angles,opt);
rec = sirtden(bb,W{1},opt.levels(1),opt.iter,opt.minc,opt.maxc);


for i =2:max(size(opt.levels))
    numgrid = size(v1{i},1);
    rec = imresize(rec,[opt.levels(i),opt.levels(i)]);
    
    for k = 1:numgrid
        for j = 1:numgrid
            
            pixels = get_pixels(v1{i}(j,:),v1{i}(k,:),opt.levels(i));
            
            temp = rec;
            temp(pixels)=0;
            b2 = bb - ...
                W{i}*reshape(temp,opt.levels(i)^2,slices);
            clear temp;
            
            temp = sirtden(b2,W{i}(:,pixels),...
                v1{i}(j,2)-v1{i}(j,1)+1,opt.iter,opt.minc,opt.maxc);
            
            rec(v2{i}(j):v2{i}(j+1)-1,v2{i}(k):v2{i}(k+1)-1,:)=...
                temp(v2{i}(j)-v1{i}(j,1)+1:v2{i}(j+1)-v1{i}(j,1),...
                v2{i}(k)-v1{i}(k,1)+1:v2{i}(k+1)-v1{i}(k),:);
        end
    end
end
            
            
	
    



end

