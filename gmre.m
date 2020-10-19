%GMRES algorithm
function x=gmre(A,b)
e=0.0001;
n=length(b);x0=zeros(n,1);r=b-A*x0;
beta=norm(r);V(:,1)=r/beta;
for m=1:n
    %do arnoldi process
    [V,Hhat]=arno(A,r,m);
    ym=lsp(Hhat,beta);
    x=x0+V(:,1:m)*ym;
    if norm(b-A*x)<e
        break;
    end
end
end