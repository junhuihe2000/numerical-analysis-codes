%conjugate gradient algorithm
function x=congra(A,b)
n=length(b);e=0.001;
x=zeros(n,1);r1=b-A*x;p=r1;
for k=0:n
    alpha=inp(r1,r1)/inp(A*p,p);
    x=x+alpha*p;
    if  inp(alpha*p,alpha*p)<e^2
        return;
    end
    r2=r1-alpha*A*p;
    beta=inp(r2,r2)/inp(r1,r1);
    p=r2+beta*p;
    if r2==0|p==0|k==n-2
        return;
    end
    r1=r2;
end
alpha=inp(r2,r2)/inp(A*p,p);
x=x+alpha*p;
end

%calculate inner product
function x=inp(a,b)
x=sum(a.*b);
end