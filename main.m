clc;
clear all;
%% input data
dt=case5;
node=length(dt.bus(:,1));
ng=length(dt.gen(:,1));
line=length(dt.branch(:,1));
Bus = dt.bus; 
Gen= dt.gen;
Gencost=dt.gencost;
Branch=dt.branch;
Vmax=max(Bus(:,12));
Vmin=min(Bus(:,13));
%%
ub = [Vmax*ones(1,node), Vmax*ones(1,node), Gen(:,9)', Gen(:,4)'];
lb = [-Vmax*ones(1,node), -Vmax*ones(1,node), Gen(:,10)', Gen(:,5)'] ;
x0=  [1.0*ones(1,node), 0.0*ones(1,node), 0.0*ones(1,ng), 0.0*ones(1,ng)]; 

A= [];
b= [];
Aeq =[];
beq =[];
options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',50000);
[x,fval,exitflag,output,lambda]=fmincon(@(x) funcname(x,node,ng,Gencost),x0,A,b,Aeq,beq,lb,ub,@(x) nonconst(x,node,ng,line,Branch,Bus,Gen),options)
x_cen = x
lam_cen = lambda.eqnonlin
fval
%%
lam = 300.0*ones(1,node);
es = 1.0*ones(1,node);
fs = 0.0*ones(1,node);
mu12=0;
exitflag=1;
output_tot =[];
for i=1:node
    output_tot = [output_tot, exitflag, es(i), fs(i), lam(i)];
end
output_tot = [output_tot;output_tot];
tic;
err=10*ones(1,node);
% for iter=3:100
iter=3;
while max(max(abs(err)))>0.001
    iter
    %%nodal optimization
    for n=1:node
        if ismember(n,Gen(:,1))
%             ub_n = [2.0, 2.0, Gen(find(Gen(:,1)==n),9), Gen(find(Gen(:,1)==n),4)];
%             lb_n = [0.0, -2.0, Gen(find(Gen(:,1)==n),10), Gen(find(Gen(:,1)==n),5)];
%             x0_n=  [1.0, 0.0, 0.0, 0.0]; 
            ub_n = [2.0, 2.0, Gen(find(Gen(:,1)==n),9), Gen(find(Gen(:,1)==n),4), Inf, Inf];
            lb_n = [0.0, -2.0, Gen(find(Gen(:,1)==n),10), Gen(find(Gen(:,1)==n),5), 0.0, -Inf];
            x0_n=  [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]; 

            A= [];
            b= [];
            Aeq =[];
            beq =[];
            options = optimoptions('fmincon','Display','off');
            [x_n,fval,exitflag_n,output,lambda_n]=fmincon(@(x) funcname_n(x,n,lam,es,fs,line,Gencost,Branch,Gen,mu12),x0_n,A,b,Aeq,beq,lb_n,ub_n,@(x) nonconst_n(x,n,es,fs,line,Branch,Bus,Gen),options);
            es(n)=x_n(1);
            fs(n)=x_n(2);
            Pg(find(Gen(:,1)==n))=x_n(3);
            Qg(find(Gen(:,1)==n))=x_n(4);
            extra(n,1)=x_n(5);
            extra(n,2)=x_n(6);
            lam(n)= lambda_n.eqnonlin(1);
%             if n==1
%                 mu12 = lambda_n.ineqnonlin(end);
%             end
            if exitflag_n~=1
                display(exitflag_n);
            end 
%             if abs(lam(n))-abs(output_tot(iter-1,n*4))>100
%                 es(n) = output_tot(iter-1,4*n-2) - (output_tot(iter-2,4*n-2)-output_tot(iter-1,4*n-2));
%                 fs(n) = output_tot(iter-1,4*n-1) - (output_tot(iter-2,4*n-1)-output_tot(iter-1,4*n-1));
%                 lam(n) = output_tot(iter-1,4*n) - (output_tot(iter-2,4*n)-output_tot(iter-1,4*n));
%             end
        else 
            ub_n = [2.0, 2.0, Inf, Inf];
            lb_n = [0.0, -2.0, 0.0, -Inf];
            x0_n=  [1.0, 0.0, 0.0, 0.0]; 
%             ub_n = [2.0, 2.0];
%             lb_n = [0.0, -2.0];
%             x0_n=  [1.0, 0.0]; 
            A= [];
            b= [];
            Aeq =[];
            beq =[];
            options = optimoptions('fmincon','Display','off');
            [x_n,fval,exitflag_n,output,lambda_n]=fmincon(@(x) funcname_n(x,n,lam,es,fs,line,Gencost,Branch,Gen,mu12),x0_n,A,b,Aeq,beq,lb_n,ub_n,@(x) nonconst_n(x,n,es,fs,line,Branch,Bus,Gen),options);
            es(n)=x_n(1);
            fs(n)=x_n(2);
            extra(n,1)=x_n(3);
            extra(n,2)=x_n(4);
            lam(n)= lambda_n.eqnonlin(1);
%             if n==1
%                 mu12 = lambda_n.ineqnonlin(end);
%             end
            if exitflag_n~=1
                display(exitflag_n);
            end 
%             if abs(lam(n))-abs(output_tot(iter-1,n*4))>100
%                 es(n) = output_tot(iter-1,4*n-2) - (output_tot(iter-2,4*n-2)-output_tot(iter-1,4*n-2));
%                 fs(n) = output_tot(iter-1,4*n-1) - (output_tot(iter-2,4*n-1)-output_tot(iter-1,4*n-1));
%                 lam(n) = output_tot(iter-1,4*n) - (output_tot(iter-2,4*n)-output_tot(iter-1,4*n));
%             end
        end 
        output_tot(iter,4*n-3) = exitflag_n;
        output_tot(iter,4*n-2) = es(n);
        output_tot(iter,4*n-1) = fs(n);
        output_tot(iter,4*n) = lam(n);
        err(n)=output_tot(iter,4*n)-output_tot(iter-1,4*n)
    end
    Poutput(iter,:)=Pg;
    Qoutput(iter,:)=Qg;
    iter=iter+1;
end
 toc;
 %voltage
 table_volt= [x_cen(1:node);
              x_cen(node+1:2*node)];
for i=1:node
    table_dvolt(1,i)=output_tot(end,4*i-2);
    table_dvolt(2,i)=output_tot(end,4*i-1);
end
% price
for i=1:node
    table_lam(1,i) = lam_cen(2*i-1);
end
for i=1:node
    table_dlam(1,i)=output_tot(end,4*i);
end
%power
table_power = [x_cen(2*node+1:2*node+2*ng)];
table_dpower = [Pg, Qg];

%%lambda plot
% figure(1)
% for i=1:14
% plot(1:length(output_tot(:,4*i)),output_tot(:,4*i))
% hold on
% end
% hold on
% axes('Position',[.7 .7 .2 .2])
% box on
% for i=1:14
%     y2=output_tot(60:end,4*i);
%     plot(60:60-1+length(y2(:,1)),y2)
%     hold on
% end


%%% Pg plot
% for i=1:ng
% plot(1:length(Poutput(:,1)),Poutput(:,i))
% hold on    
% end

%%% Qg plot
% for i=1:ng
% plot(1:length(Qoutput(:,1)),Qoutput(:,i))
% hold on    
% end

%%% volt plot
% for i=1:14
% tx1(i)=output_tot(end,4*i-2);
% tx2(i)=output_tot(end,4*i-1);
% end
% figure
% scatter(tx1,tx2)
% hold on
% syms re1 im1
% ezplot(re1^2 + im1^2 ==1.21,[0.8,1.1,-0.4,0.4])
% hold on
% ezplot(re1^2 + im1^2 ==0.81,[0.8,1.1,-0.4,0.4])
% 


