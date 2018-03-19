function f = funcname(x,node,ng,Gencost)
f=0;
for i=1:ng
    f = f + Gencost(i,5)*x(2*node+i)*x(2*node+i)*100 + Gencost(i,6)*x(2*node+i)*10;
% f = f + Gencost(i,6)*x(2*node+i);
end
end