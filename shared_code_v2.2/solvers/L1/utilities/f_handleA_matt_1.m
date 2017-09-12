function y = f_handleA_matt_1(u,mode,st,N)


switch mode
    case 1
        F = reshape(u,2*N+1,2*N+1);
        y = nufft(F.',st);
    case 2
        y = nufft_adj(u,st).';
	y = y(:);
    otherwise
        error('Unknown mode passed to f_handleA_matt_1 in ftv_cs.m');
end

end
