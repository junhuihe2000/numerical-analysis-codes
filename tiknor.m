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
