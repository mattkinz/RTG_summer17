function y = combo_operator(A1,A1t,A2,A2t,x,mode,A1dim)

switch mode
    case 1
        y1 = A1(x);
        y2 = A2(x);
        y = [y1;y2];
    case 2
        y1 = A1t(x(1:A1dim));
        y2 = A2t(x(A1dim+1:end));
        y = [y1;y2];        
end

end