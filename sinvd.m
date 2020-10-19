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


