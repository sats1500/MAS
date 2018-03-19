function [c,ceq] = nonconst_n(x,n,es,fs,line,Branch,Bus,Gen)
Branch(:,14)=complex(Branch(:,3),Branch(:,4));
Branch(:,15)=1./Branch(:,14);
P12max=0.50;
ceq=[];
c=[];

ceq(1,1)= Bus(n,3);
ceq(2,1)= Bus(n,4);

if ismember(n,Gen(:,1))
    ceq(1,1)= ceq(1,1) - x(3);
    ceq(2,1)= ceq(2,1) - x(4);
end

for j=1:line
    n1=Branch(j,1);
    n2=Branch(j,2);
    if n==n1 
        ceq(1,1)= ceq(1,1) + (x(1)^2 + x(2)^2 - x(1)*es(n2) - x(2)*fs(n2))*real(Branch(j,15)) - (x(1)*fs(n2) - x(2)*es(n2))*imag(Branch(j,15)) ;
        ceq(2,1)= ceq(2,1) + (x(1)^2 + x(2)^2 - x(1)*es(n2) - x(2)*fs(n2))*imag(Branch(j,15)) + (x(1)*fs(n2) - x(2)*es(n2))*real(Branch(j,15)) ;
    end
    if n==n2
        ceq(1,1)= ceq(1,1) + (x(1)^2 + x(2)^2 - x(1)*es(n1) - x(2)*fs(n1))*real(Branch(j,15)) - (x(1)*fs(n1) - x(2)*es(n1))*imag(Branch(j,15)) ;
        ceq(2,1)= ceq(2,1) + (x(1)^2 + x(2)^2 - x(1)*es(n1) - x(2)*fs(n1))*imag(Branch(j,15)) + (x(1)*fs(n1) - x(2)*es(n1))*real(Branch(j,15)) ;
    end
       
end

% for j=1:line
%     n1=Branch(j,1);
%     n2=Branch(j,2);
%     if n==n1 
%         ceq(1,1)= ceq(1,1) - (es(n2)^2 + fs(n2)^2 - es(n2)*x(1) - fs(n2)*x(2))*real(Branch(j,15)) + (es(n2)*x(2) - fs(n2)*x(1))*imag(Branch(j,15)) ;
%         ceq(2,1)= ceq(2,1) - (es(n2)^2 + fs(n2)^2 - es(n2)*x(1) - fs(n2)*x(2))*imag(Branch(j,15)) - (es(n2)*x(2) - fs(n2)*x(1))*real(Branch(j,15)) ;
%     end
%     if n==n2
%         ceq(1,1)= ceq(1,1) - (es(n1)^2 + fs(n1)^2 - es(n1)*x(1) - fs(n1)*x(2))*real(Branch(j,15)) + (es(n1)*x(2) - fs(n1)*x(1))*imag(Branch(j,15)) ;
%         ceq(2,1)= ceq(2,1) - (es(n1)^2 + fs(n1)^2 - es(n1)*x(1) - fs(n1)*x(2))*imag(Branch(j,15)) - (es(n1)*x(2) - fs(n1)*x(1))*real(Branch(j,15)) ;
%     end
%        
% end
ceq(1,1)=ceq(1,1)-x(end-1);
ceq(2,1)=ceq(2,1)-x(end);


c(1,1)=x(1)^2 + x(2)^2 - Bus(n,12)^2;
c(2,1)= -x(1)^2 - x(2)^2 + Bus(n,13)^2;
% if n==1
%     c(end+1,1)=(x(1)^2 + x(2)^2 - x(1)*es(2) - x(2)*fs(2))*real(Branch(1,15)) - (x(1)*fs(2) - x(2)*es(2))*imag(Branch(1,15)) - P12max ;
% end




if Bus(n,2)==3
    ceq(end+1,1) = x(2);
    ceq(end+1,1) = x(1)-1.0;
end
    
end