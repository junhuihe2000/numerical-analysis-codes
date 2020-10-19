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

function x=inp(a,b)
x=sum(a.*b);
end