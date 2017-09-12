function x = RW_LS(x,mode,A,R)

switch mode
    case 1
        x = A(x,1);
        x = x.*R;
    case 2     
        x = x.*R;
        x = A(x,2);
end