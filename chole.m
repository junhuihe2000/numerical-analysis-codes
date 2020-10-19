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
