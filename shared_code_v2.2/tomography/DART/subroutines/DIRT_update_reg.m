function Uchunk = DIRT_update_reg(W,M,bchunk,Uchunk,set,fixed,update_type,opt,n)

        x = cell(size(Uchunk,2),1);
      
        
        for j = 1:size(Uchunk,2)
            bchunk(:,j) = bchunk(:,j) - W(:,fixed{j})*Uchunk(fixed{j},j);
        end
        
        if strcmp(update_type,'CGLS')
            
            for j = 1:size(Uchunk,2)
                Wt=W;
                Wt(:,fixed{j})=0;
                x{j} = cgls([W;M],[bchunk(:,j);zeros(size(M,1),1)],0,10^(-6),opt.inner_iter);
            end
        else
            x = sirt_sub_chunk_reg(W,M,bchunk,opt.inner_iter,Uchunk,...
                set,fixed,opt.minc,opt.grays,opt.lambda1,opt.lambda2);
        end
        
        for j = 1:size(Uchunk,2)
            Uchunk(:,j)=x(:,j);
        end
        Uchunk = reshape(Uchunk,n,n,size(Uchunk,2));
        
end