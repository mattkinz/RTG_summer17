function Uchunk = DIRT_update(W,Wn,bchunk,Uchunk,set,fixed,update_type,opt,c,n)

        x = cell(size(Uchunk,2),1);
        inits = cell(size(Uchunk,2),1);
        cs = cell(size(Uchunk,2),1);
        
        for j = 1:size(Uchunk,2)
            bchunk(:,j) = bchunk(:,j) - W(:,fixed{j})*Uchunk(fixed{j},j);
            if strcmp(update_type,'SIRT');
                inits{j}=Uchunk(set{j},j);
                %cs{j}=c(set{j});
            end
        end
        
        if strcmp(update_type,'CGLS')
            for j = 1:size(Uchunk,2)
                x{j} = cgls(W(:,set{j}),bchunk(:,j),0,10^(-6),opt.inner_iter);
            end
        else
            x = sirt_sub_chunk(W,Wn,bchunk,opt.inner_iter,inits,set,c,opt.minc,opt.grays);
        end
        
        for j = 1:size(Uchunk,2)
            Uchunk(set{j},j)=x{j};
        end
        Uchunk = reshape(Uchunk,n,n,size(Uchunk,2));
        
end