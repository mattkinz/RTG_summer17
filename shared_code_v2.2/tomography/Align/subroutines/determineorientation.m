function [stack,shifts,rotation] = determineorientation(stack,shifts,opt)

%subroutine for autoalign


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

[m,n,k]=size(stack);
st = zeros(n,m,k);
for i =1:k
    st(:,:,i)=imrotate(stack(:,:,i),90);
end
fprintf('Weighting the projections\n');
stw = weighting2(st,opt.window_accuracy);
stackw = weighting2(stack,opt.window_accuracy);


fprintf('Determining if the stack should be rotated.\n');
fprintf('(note: the tilt axis should be horizontal)\n');

s1 = sum(stackw,1);
s2 = sum(stw,1);


mu1 = median(s1,3);
mu2 = median(s2,3);
b1 = zeros(n,1);
b2 = zeros(m,1);

for i =1:n
    b1(i)=norm(reshape(s1(1,i,:)-mu1(i),size(s1,3),1),1);
end
for i = 1:m
    b2(i)=norm(reshape(s2(1,i,:)-mu2(i),size(s2,3),1),1);
end

b1 = norm(b1,1);
b2 = norm(b2,1);
[~,which]=min([b1,b2]);
if which==2
    stack = st;
    stackw = stw;
    n=size(stack,2);
    shiftt=shifts;
    shifts(:,1)=-shiftt(:,2);
    shifts(:,2)=shiftt(:,1);
    fprintf('Determined that stack needed to be rotated\n\n');
    rotation = true;
else
    fprintf('No rotation needed\n\n');
    rotation=false;
end

clear stw;clear st;clear shiftt;
close(figure(1));

if strcmp(opt.xalign_type,'local')
    fprintf('Performing local x-direction alignment, \n using conservation of mass\n\n');
    [stack,v,stackw] = ...
        corrx_local(stack,[round(n/4),round(3*n/4)],stackw);
    shifts(:,2)=shifts(:,2)+v;

    for i = 1:5    
        if max(v)~=20
            break;
        else
            fprintf('Reiterating local x-direction alignment\n');
            [stack,v,stackw] = ...
        corrx_local(stack,[round(n/4),round(3*n/4)],stackw);
        shifts(:,2)=shifts(:,2)+v;
        end
    end
else
    fprintf('Performing global x-direction alignment, \n using conservation of mass\n');
    [stack,v,stackw] = ...
        corrx_global(stack,[round(n/4),round(3*n/4)],stackw);
    shifts(:,2)=shifts(:,2)+v;

    for i = 1:5    
        if max(v)<10
            break;
        else
            fprintf('Reiterating global x-direction alignment\n');
            [stack,v,stackw] = ...
        corrx_global(stack,[round(n/4),round(3*n/4)],stackw);
        shifts(:,2)=shifts(:,2)+v;
        end
    end
end
end




