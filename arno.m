%Arnoldi process
function [V,Hhat]=arno(A,r,m)
n=length(A);V=zeros(n,m+1);H=zeros(m+1);
V(:,1)=r/norm(r);
for k=1:m
    t=A*V(:,k);
    for i=1:k
        H(i,k)=inp(A*V(:,k),V(:,i));
        t=t-H(i,k)*V(:,k);
    end
    H(k+1,k)=norm(t);
    V(:,k+1)=t/norm(t);
end
em=zeros(1,m);em(m)=1;
Hhat=[H(1:m,1:m);H(m+1,m)*em];
end

function x=inp(a,b)
x=sum(a.*b);
end
    