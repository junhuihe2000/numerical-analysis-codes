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