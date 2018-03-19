function [c,ceq] = nonconst(x,node,ng,line,Branch,Bus,Gen)
Branch(:,14)=complex(Branch(:,3),Branch(:,4));
Branch(:,15)=1./Branch(:,14);
P12max=0.50;
ceq=[];
c=[];
for i=1:node
    ceq(2*i-1,1)= Bus(i,3);
    ceq(2*i,1)= Bus(i,4);
    if ismember(i,Gen(:,1))
        ceq(2*i-1,1)= ceq(2*i-1,1)-x(2*node + find(Gen(:,1)==i));
        ceq(2*i,1)= ceq(2*i,1)-x(2*node + ng + find(Gen(:,1)==i));
    end
    for j=1:line
        if i==Branch(j,1) 
            ceq(2*i-1,1)= ceq(2*i-1,1) + (x(i)^2 + x(node+i)^2 - x(i)*x(Branch(j,2)) - x(node+i)*x(node+Branch(j,2)))*real(Branch(j,15)) - (x(i)*x(node+Branch(j,2)) - x(node+i)*x(Branch(j,2)))*imag(Branch(j,15)) ;
            ceq(2*i,1)= ceq(2*i,1) + (x(i)^2 + x(node+i)^2 - x(i)*x(Branch(j,2)) - x(node+i)*x(node+Branch(j,2)))*imag(Branch(j,15)) + (x(i)*x(node+Branch(j,2)) - x(node+i)*x(Branch(j,2)))*real(Branch(j,15)) ;
        end
        if i==Branch(j,2)
            ceq(2*i-1,1)= ceq(2*i-1,1) + (x(i)^2 + x(node+i)^2 - x(i)*x(Branch(j,1)) - x(node+i)*x(node+Branch(j,1)))*real(Branch(j,15)) - (x(i)*x(node+Branch(j,1)) - x(node+i)*x(Branch(j,1)))*imag(Branch(j,15)) ;
            ceq(2*i,1)= ceq(2*i,1) + (x(i)^2 + x(node+i)^2 - x(i)*x(Branch(j,1)) - x(node+i)*x(node+Branch(j,1)))*imag(Branch(j,15)) + (x(i)*x(node+Branch(j,1)) - x(node+i)*x(Branch(j,1)))*real(Branch(j,15)) ;
        end
    end
    c(2*i-1,1)=x(i)^2 + x(node+i)^2 - Bus(i,12)^2;
    c(2*i,1)= -x(i)^2 - x(node+i)^2 + Bus(i,13)^2;
    
end
% c(end+1,1)=(x(1)^2 + x(node+1)^2 - x(1)*x(2) - x(node+1)*x(node+2))*real(Branch(1,15)) - (x(1)*x(node+2) - x(node+1)*x(2))*imag(Branch(1,15)) - P12max ;


ceq(2*node+1)=x(node+Bus(find(Bus(:,2)==3),1));
ceq(2*node+2)=x(Bus(find(Bus(:,2)==3),1))-1.0;

% temp1= x(Branch(:,1)).^2 ;
% temp2= x(node+Branch(:,1)).^2 ;
% temp3= x(Branch(:,1)).*x(Branch(:,2)) ;
% temp4= x(node+Branch(:,1)).*x(node+Branch(:,2));
% temp5= real(Branch(:,15));
% temp6= x(Branch(:,1)).*x(node+Branch(:,2)) ;
% temp7= x(node+Branch(:,1)).*x(Branch(:,2));
% temp8= imag(Branch(:,15));
% 
% temp9= x(Branch(:,2)).^2 ;
% temp10= x(node+Branch(:,2)).^2 ;
% temp11= x(Branch(:,2)).*x(Branch(:,1)) ;
% temp12= x(node+Branch(:,2)).*x(node+Branch(:,1));
% temp13= real(Branch(:,15)) ;
% temp14= x(Branch(:,2)).*x(node+Branch(:,1));
% temp15= x(node+Branch(:,2)).*x(Branch(:,1));
% temp16= imag(Branch(:,15));
% Ppowerflow12=sum((temp1(:) + temp2(:) - temp3(:) - temp4(:)).*temp5 - (temp6(:) - temp7(:)).*temp8) ;
% Ppowerflow21=sum((temp9(:) + temp10(:) - temp11(:) - temp12(:)).*temp13 - (temp14(:) - temp15(:)).*temp16);
% Qpowerflow12=sum((temp1(:) + temp2(:) - temp3(:) - temp4(:)).*temp8 + (temp6(:) - temp7(:)).*temp5) ;
% Qpowerflow21=sum((temp9(:) + temp10(:) - temp11(:) - temp12(:)).*temp16 + (temp14(:) - temp15(:)).*temp13) ;
% ceq(2*node+3)=sum(Bus(:,3))-sum(x(2*node:2*node+ng)) + Ppowerflow12 + Ppowerflow21 ;
% ceq(2*node+4)=sum(Bus(:,4))-sum(x(2*node:2*node+2*ng)) + Qpowerflow12 + Qpowerflow21 ;

end