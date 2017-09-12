function y = sub_oper(W,x,mode,set)


%Written by: Toby Sanders
%Comp. & Applied Math Dept., Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016

slices = numel(set);
[m,n]=size(W);


switch mode
    case 1
        x = reshape(x,n,slices);
        y = zeros(m,slices);
        for i = 1:slices
            y(:,i)=W(:,set{i})*x(set{i},i);
        end
        y = y(:);
    case 2
        y = zeros(n,slices);
        c=0;
        x = reshape(x,m,slices);
        for i = 1:slices
            y(set{i},i)=...
                (x(:,i)'*W(:,set{i}))';
        end
        y = y(:);
        
end