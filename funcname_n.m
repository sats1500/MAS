function f = funcname_n(x,n,lam,es,fs,line,Gencost,Branch,Gen,mu12)
Branch(:,14)=complex(Branch(:,3),Branch(:,4));
Branch(:,15)=1./Branch(:,14);
f=0;
if ismember(n,Gen(:,1))
    f = f + Gencost(find(Gen(:,1)==n),5)*x(3)*x(3)*100 + Gencost(find(Gen(:,1)==n),6)*x(3)*10;
%     f = f + Gencost(find(Gen(:,1)==n),5)*x(3) ;
end 
mu12=0;

for j=1:line
    n1=Branch(j,1);
    n2=Branch(j,2);
    if n==n1
        f = f + lam(n2)*((es(n2)^2 + fs(n2)^2 - es(n2)*x(1) - fs(n2)*x(2))*real(Branch(j,15)) - (es(n2)*x(2) - fs(n2)*x(1))*imag(Branch(j,15)));
    end 
    if n==n2
        if n==2
            f = f + (lam(n1)+mu12)*((es(n1)^2 + fs(n1)^2 - es(n1)*x(1) - fs(n1)*x(2))*real(Branch(j,15)) - (es(n1)*x(2) - fs(n1)*x(1))*imag(Branch(j,15)));
        else
            f = f + lam(n1)*((es(n1)^2 + fs(n1)^2 - es(n1)*x(1) - fs(n1)*x(2))*real(Branch(j,15)) - (es(n1)*x(2) - fs(n1)*x(1))*imag(Branch(j,15)));
        end
        
    end 
end 

f = f + 2*max(abs(lam))*x(end-1) + 2*max(abs(lam))*x(end)*x(end);
% for j=1:line
%     n1=Branch(j,1);
%     n2=Branch(j,2);
%     if n==n1
%         f = f - lam(n2)*((x(1)^2 + x(2)^2 - es(n2)*x(1) - fs(n2)*x(2))*real(Branch(j,15)) - (fs(n2)*x(1) - es(n2)*x(2))*imag(Branch(j,15)));
%     end 
%     if n==n2
%         if n==2
%             f = f - (lam(n1)+mu12)*((x(1)^2 + x(2)^2 - es(n1)*x(1) - fs(n1)*x(2))*real(Branch(j,15)) - (fs(n1)*x(1) - es(n1)*x(2))*imag(Branch(j,15)));
%         else
%             f = f - lam(n1)*((x(1)^2 + x(2)^2 - es(n1)*x(1) - fs(n1)*x(2))*real(Branch(j,15)) - (fs(n1)*x(1) - es(n1)*x(2))*imag(Branch(j,15)));
%         end
%         
%     end 
% end 


end
