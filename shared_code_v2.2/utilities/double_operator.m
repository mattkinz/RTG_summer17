function y = double_operator(A,B,u,mode)

switch mode
    case 1
        y = B(u,1);
        y = A(y,1);
    case 2
        y = A(u,2);
        y = B(y,2);
end

        