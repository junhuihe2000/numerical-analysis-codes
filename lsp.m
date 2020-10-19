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