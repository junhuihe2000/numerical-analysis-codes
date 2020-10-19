%solve a problem about linear equation
%input a Hilbert matrix
n=15;
H=Hilm(n);
x=ones(n,1);
b=H*x;
%gauss elimination
y1=gaueli(H,b)
%cholesky decompose
y2=chole(H,b)
%tikhonov normalized
y3=tiknor(H,b)
%conjugate gradient
y4=congra(H,b)
%GMRES
y5=gmre(H,b)

%Gauss column pivot elimination
function y=gaueli(A,b)
n=length(A);y=zeros(n,1);B=[A,b];
%elimination
for k=1:n-1
    m=k;
    for l=k:n
        if abs(B(l,k))>abs(B(m,k))
            m=l;
        end
    end
    B([k m],:)=B([m k],:);
    for i=k+1:n
        l=B(i,k)/B(k,k);
        B(i,:)=B(i,:)-B(k,:)*l;
    end
end
%back substitution
A=B(:,1:n);b=B(:,n+1);
y(n)=b(n)/A(n,n);
for i=n-1:-1:1
    y(i)=b(i);
    for j=n:-1:i+1
        y(i)=y(i)-A(i,j)*b(j);
    end
    y(i)=y(i)/A(i,i);
end
end

%Cholesky decompose to solve linear qauation
function x=chole(A,b)
n=length(A);y=zeros(n,1);x=y;
%Cholesky decompose
d=zeros(n,1);T=zeros(n);L=T;
d(1)=A(1,1);
for i=2:n
    T(i,1)=A(i,1);L(i,1)=T(i,1)/d(1);m=0;
    for j=2:i-1
        for k=1:j-1
            m=m+T(i,k)*L(j,k);
        end
        T(i,j)=A(i,j)-m;
        L(i,j)=T(i,j)/d(j);
    end
    m=0;
    for k=1:i-1
        m=m+T(i,k)*L(i,k);
    end
    d(i)=A(i,i)-m;
end
%solve the equation Ly=b
y(1)=b(1);
for i=2:n
    m=0;
    for k=1:i-1
        m=m+L(i,k)*y(k);
    end
    y(i)=b(i)-m;
end
%solve the equation DL*x=y
x(n)=y(n)/d(n);
for i=n-1:-1:1
    m=0;
    for k=i+1:n
        m=m+L(k,i)*x(k);
    end
    x(i)=y(i)/d(i)-m;
end
end

%Tikhonov normalized algorithm
function x=tiknor(A,b)
[U,D,V]=sinvd(A);d=diag(D);n=length(b);
r=rank(D);
%a is the normalized constant,1>>a>0
a=min(100*d(1)*d(r),0.0001);
x=zeros(n,1);
for j=1:r
    f=d(j)/(a+d(j)^2);
    g=sum(b.*V(:,j));
    t=f*g;
    x=x+t*U(:,j);
    %x=x+(d(j)/(a+d(j)^2))*(sum(b.*V(:,j)))*U(:,j);
end
end

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

%GMRES algorithm
function x=gmre(A,b)
e=0.001;
n=length(b);x0=zeros(n,1);r=b-A*x0;
beta=norm(r);V(:,1)=r/beta;
for m=1:n
    %do arnoldi process
    [V,Hhat]=arno(A,r,m);
    ym=lsp(Hhat,beta);
    x=x0+V(:,1:m)*ym;
    if norm(b-A*x)<e
        return;
    end
end
end

%generate a Hilbert matrix
function H=Hilm(n)
H=zeros(n);
for i=1:n
    for j=1:n
        H(i,j)=1/(i+j-1);
    end
end
end

%singular value decompose
function [U,D,V]=sinvd(A)
%A=VDU'
n=length(A);
M=A'*A;
[Us,Ds]=eig(M);
[d,ind]=sort(diag(Ds),'descend');
Ds=Ds(ind,ind);Us=Us(:,ind);
U=gsor(Us);
N=A*A';
[Vs,Ds]=eig(N);
[d,ind]=sort(diag(Ds),'descend');
Ds=Ds(ind,ind);Vs=Vs(:,ind);
Ds=abs(Ds);
V=gsor(Vs);D=sqrtm(Ds);
d=sqrt(d);
end

%solve the least square problem
%min||beta*e1-H*x||
function x=lsp(H,beta)
n=length(H);e1=zeros(n,1);e1(1)=1;
%solve the equation H'Hy=beta*H'e1
M=H'*H;b=beta*H'*e1;
if n<=2
    x=b/M;
    return;
end
x=congra(M,b);
end

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

%Gram-Schmidt normalized
function A=gsor(B)
n=length(B);A=B;
A(:,1)=A(:,1)/norm(A(:,1));
for k=2:n
    for i=1:k-1
        A(:,k)=A(:,k)-inp(A(:,k),A(:,i))*A(:,i);
    end
    A(:,k)=A(:,k)/norm(A(:,k));
end
end