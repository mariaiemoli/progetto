%prova di intersezione di due fratture 1D

close all
clear all
clc
format

lista_fratture_param4

%calcolo delle intersezioni
fx=zeros(Nf);
progressivo=zeros(Nf);
t_int=zeros(Nf);
Nj=zeros(Nf);
for i=1:Nf
    quali{i}=[];
   for j=1:Nf
       if (i~=j)
    %   valuto le due sugli stessi nodi e cerco in quale elemento DELLA PRIMA si incrociano
%        dx=(x_map{i}(tt{i})-x_map{j}(tt{i}));
%        [blabla,Nj(i,j)]=min(dx(1,:).^2+dx(2,:).^2);
%        
%       xx=[t{i}(Nj(i,j)):(t{i}(Nj(i,j)+1)-t{i}(Nj(i,j)))/100:t{i}(Nj(i,j)+1)];  %%raffino solo lì e trovo il punto esatto
%        if (xx(1)<0 || xx(1)<0 || xx(end)>1 ||xx(end)>1)
%            Nj(i,j)=0;
%        else
%        dx=(x_map{i}(xx)-x_map{j}(xx));
%        [blabla,jj]=min(dx(1,:).^2+dx(2,:).^2);
%        t_int(i,j)=xx(jj);
%    
%         %cerco l'intersezione con newton
        toll=1.0e-6;
        err=toll+1;
        maxit =100;
        t_sol=[0.5;0.5];
        nit=0;
       while err>toll && nit<maxit
           nit=nit+1;
           res=x_map{i}(t_sol(1))-x_map{j}(t_sol(2));
           jac=[Tau{i}(t_sol(1)), -Tau{j}(t_sol(2))];
           if (det(jac)==0) break; 
           end
           dx=-jac\res;
           t_sol=t_sol+dx;
           err=norm(res);
       end
       if (err>toll )
        toll=1.0e-6;
        err=toll+1;
        maxit =100;
        t_sol=[1;1];
        nit=0;
       while err>toll && nit<maxit
           nit=nit+1;
           res=x_map{i}(t_sol(1))-x_map{j}(t_sol(2));
           jac=[Tau{i}(t_sol(1)), -Tau{j}(t_sol(2))];
           if (det(jac)==0) break;
           end
           dx=-jac\res;
           t_sol=t_sol+dx;
           err=norm(res);
       end
       end
       if (err>toll)
       t_sol=-1;
       end
       if (t_sol(1)<0 || t_sol(1)>1)
           N_j(i,j)=0;
       else
       t_int(i,j)=t_sol(1);
       Nj(i,j)=ceil(t_sol(1)/h(i));
       fx(i,j)=(t_int(i,j)-h(i)*(Nj(i,j)-1))/h(i);   %espresso come frazione dell'elemento
       quali{i}=[quali{i} j];
       KI{i,j}=(power(i)==power(j))*(eta(i) + eta(j))^-1*eye(2) + (power(i)>power(j))*eta(i)^-1*eye(2)+  (power(i)<power(j))*eta(j)^-1*eye(2);
%       KI{i,j}=10^7*eye(2);
       I(i,j)=d(i)*d(j)/(1-(tau{i}(t_int(i,j))'*tau{j}(t_int(i,j)))^2)^0.5;
       end
       end
 
   end
end

%conto il numero di intersezioni

Nint=sum(Nj>0,2);
Nint_tot=sum(Nint);

%decido una numerazione globale delle intersezioni 
global_id=zeros(Nf);
count=1;
for i=1:Nf
    for j=i+1:Nf
        if (Nj(i,j)>0 )
            global_id(i,j)=count;
            count=count+1;
        end
    end
end

%faccio in modo che sia "simmetrica" cioè se sono 4 nodi di intersezione
% p1 corrisponde a p5 e così via
for i=1:Nf
    for j=1:i-1
        if (Nj(i,j)>0 )
            global_id(i,j)=global_id(j,i) + Nint_tot/2;
        end
    end
end


%calcolo il numero di dof per comodità nell'assemblaggio

dof=2*N + ones(size(N)) + 3*Nint';
dofP=N + Nint';
dofV=N + ones(size(N)) + 2*Nint';

dof_tot = sum(dof);

M=zeros(dof_tot+Nint_tot,dof_tot+Nint_tot);
F=zeros(dof_tot+Nint_tot,1);

%assemblo matrici e termine noto per le fratture. 
Tshift=[0 cumsum(dof)];
shift=0;

for i=1:Nf
%assemblo
ti=tt{i};
[A,B,C,rhs,NO]=solver_darcy_1d_multi_disc_KIparam(N(i),Nj(i,:),fx(i,:),fx(:,i),i,quali{i},Nint(i),h(i).*ds{i}(ti),source{i}(ti),...
    p_bc{i},v_bc{i},t_bc{i},eta, d,KI, tau,I,t_int);
Mi=[A -B;B' zeros(N(i)+Nint(i),N(i)+Nint(i))];  %costruisco la m e la copio
M(shift+1:shift+dof(i),shift+1:shift+dof(i))=Mi;

%ora i termini di accoppiamento. per ogni frattura ne ho più di uno
for k=1:length(quali{i})
    prog=min(global_id(i,quali{i}(k)), global_id(quali{i}(k),i));  %primo e secondo indice dell'intersezione
    prog2=max(global_id(i,quali{i}(k)), global_id(quali{i}(k),i));
    M(shift+1:shift+dofV(i),end-Nint_tot+prog)=C(:,k);  %questo termine è p_hat vdotn 

Njj=(setdiff(Nj(:,quali{i}(k)),0));
[blabla]=ismember(Njj, Nj(i,quali{i}(k)));
[blabla2,where]=max(blabla);
extV1altra=N(quali{i}(k))+2+(where-1)*2;
extV2altra=N(quali{i}(k))+3+(where-1)*2;

M(shift+1:shift+dofV(i),Tshift(quali{i}(k))+Nj(quali{i}(k),i))=NO{k}(:,1);  %questo termine è p_hat vdotn 
M(shift+1:shift+dofV(i),Tshift(quali{i}(k))+Nj(quali{i}(k),i)+1)=NO{k}(:,2);  %questo termine è p_hat vdotn 
M(shift+1:shift+dofV(i),Tshift(quali{i}(k))+extV1altra)=NO{k}(:,3);  %questo termine è p_hat vdotn 
M(shift+1:shift+dofV(i),Tshift(quali{i}(k))+extV2altra)=NO{k}(:,4);  %questo termine è p_hat vdotn 


M(end-Nint_tot+prog, shift+1:shift+dofV(i))=M(end-Nint_tot+prog, shift+1:shift+dofV(i)) + C(:,k)'; %il suo simmetrico (somma flussi=0)
M(end-Nint_tot/2+prog, end-Nint_tot+prog)=1;   % e con questo forzo p1=p5, p2=p6 e così via
M(end-Nint_tot/2+prog, end-Nint_tot+prog2)=-1;

end
%termine noto
F(shift+1:shift+dof(i))=rhs;
shift=shift+dof(i);  %incremento lo shift 
end
xx=M\F;  %risolvo

figure

%per esportare la soluzione faccio un ciclo sulle fratture 
shift=0;
for i=1:Nf
   sol=xx(shift+1:shift+dof(i));
   shift=shift+dof(i);
   U{i}=sol(1:dofV(i));
   p=sol(end-dofP(i)+1:end);

   for k=1:length(quali{i})
    
         prog2=max(global_id(i,quali{i}(k)), global_id(quali{i}(k),i));
  x_int=x_map{i}(t_int(i,quali{i}(k)));

if (abs(x_int(1)-0.5)<1.0e-5 && abs(x_int(2)-0.5)<1.0e-5)
    disp('la sol è')
    xx(end-Nint_tot+prog2)
end
   plot3(x_int(1),x_int(2), xx(end-Nint_tot+prog2),'o','color',colori{i} ,'Markersize',6,'linewidth',2)
   plot3(x_int(1),x_int(2), xx(end-Nint_tot+prog2),'.','color',colori{quali{i}(k)},'Markersize',12,'linewidth',2)

   end
    hold on
   ssx=[tt{i}  tt{i}(Nj(i,quali{i}))+h(i)/100];  %metto in coda i baricentri degli elementi tagliati + una frazione micro di h per poterli ordinare
   [ssx,I]=sort(ssx);  %ordino e poi con lo stesso vettore di indici riordino la soluzione
   p=p(I);

   nodi=[t{i}  t{i}(Nj(i,quali{i}))   t{i}(Nj(i,quali{i})+1) t_int(i,quali{i})];
   nodi=unique(nodi);
   nodi=sort(nodi);
   nnodi=[nodi(1:end-1);nodi(2:end)];
   xxx=x_map{i}(nnodi);
   %plot3(xxx(1:2,:), xxx(3:4,:), (ones(2,1)*p'),'color',colori{i},'linewidth',4)
   c=colormap('jet');
   for i=1:length(p)
       k=floor(p(i)*64)+1;
       plot3(xxx(1:2,i), xxx(3:4,i), (ones(2,1)*p(i)),'color',c(k,:),'linewidth',4)
       hold on
 
   end
   
end

%calcolo K upscalata
fluxinh=0;
fluxouth=0;
fluxinv=0;
fluxoutv=0;

for i=1:Nhor
fluxinh=fluxinh+U{i}(1); 
fluxouth=fluxouth+U{i}(N(i)+1);
end
for i=1+Nhor:Nhor+Nver
fluxinv=fluxinv+U{i}(1) ;
fluxoutv=fluxoutv+U{i}(N(i)+1);
end
Kxx=fluxouth+fluxinh
Kxy=fluxoutv+fluxinv

Kyx=fluxouth+fluxinh
Kyy=fluxoutv+fluxinv




%axis([-0.5 0.5 -1 1.5 0 1])
