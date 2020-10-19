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