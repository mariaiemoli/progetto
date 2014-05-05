function [A,B,C,F,NO]=solver_darcy_1d_multi_disc(N,Nj, fx,fx_altro,chi,quali, Nint, h,source,p_b,v_bc,t_bc,eta,d, KII, tau,II,t_int);

gamma=50;
csi0=0.5;
eta_t=eta(chi)/d(chi);
%se diventa troppo lento è da cambiare: definisco la massa (tau-tau) in modo
%simbolico per integrare automaticamente
%syms xx;
%massa= [(1-xx)^2, (1-xx)*xx; xx*(1-xx),xx*xx ];

A=zeros(N+1+2*Nint,N+1+2*Nint);
F=zeros(2*N+1+3*Nint,1);    %inizializzo il termine noto

%per assemblare creo una lista degli elementi non tagliati

N_un=setdiff([1:N], Nj);

%A negli elementi non tagliati
for i=N_un
    A(i:i+1,i:i+1)=A(i:i+1,i:i+1)+[1/3 1/6;1/6 1/3]*eta_t*h(i);
end

%stessa cosa per la B
B=[zeros(N,1) eye(N,N)]-[eye(N,N) zeros(N,1)];
B=[B' zeros(N+1,Nint);zeros(2*Nint,N+Nint)];

C=zeros(N+1+2*Nint,Nint);

%termine sorgente SCALARE per il termine noto
source=[source.*h, zeros(1,Nint)];

for k=1:Nint
I=II(chi,quali(k));
KI=KII{chi,quali(k)};
kk=quali(k);   
t_chi=t_int(chi,quali(k));
t_kk=t_int(quali(k),chi);

seno=(1-(tau{chi}(t_chi)'*tau{kk}(t_kk))^2)^0.5;

quali
Ain= eta_t*h(Nj(quali(k)))*0.5*([(1)^2, 0; 0,0]+[(1-fx(kk))^2, (1-fx(kk))*fx(kk); fx(kk)*(1-fx(kk)),fx(kk)*fx(kk)]);
Aout= eta_t*h(Nj(quali(k)))*0.5*([0, 0; 0,1]+[(1-fx(kk))^2, (1-fx(kk))*fx(kk); fx(kk)*(1-fx(kk)),fx(kk)*fx(kk)]);
%Ain=eta_t*h(Nj(quali(k)))*int(massa,0 ,fx(kk));   %integro a sinistra e destra del taglio
%Aout=eta_t*h(Nj(quali(k)))*int(massa, fx(kk),1);

extV1=N+2+(k-1)*2;
extV2=N+3+(k-1)*2;
extP=N+k;

%aggiungo il contributo dell'elemento tagliato
A(Nj(kk),Nj(kk))=A(Nj(kk),Nj(kk))+Ain(1,1);
A(extV1,extV1)=A(extV1,extV1)+Aout(1,1);

A(Nj(kk),extV2)=A(Nj(kk),extV2)+Ain(1,2);
A(extV1,Nj(kk)+1)=A(extV1,Nj(kk)+1)+Aout(1,2);

A(Nj(kk)+1,extV1)=A(Nj(kk)+1,extV1)+Aout(1,2);
A(extV2,Nj(kk))=A(extV2,Nj(kk))+Ain(1,2);

A(Nj(kk)+1,Nj(kk)+1)=A(Nj(kk)+1,Nj(kk)+1)+Aout(2,2);
A(extV2,extV2)=A(extV2,extV2)+Ain(2,2);


%aggiungo alla B i contributi destro e sinistro 
B(:,Nj(kk))=0;
B(Nj(kk):Nj(kk)+1,Nj(kk))=[-fx(kk);0];
B(extV1:extV2,Nj(kk))=[0; fx(kk)];
B(Nj(kk):Nj(kk)+1,extP)=[0; 1-fx(kk)];
B(extV1:extV2,extP)=[-(1-fx(kk));0];

%valuto tau-tau nel punto di taglio per calcolare media-media e salto-salto
xx=fx(kk);

coeff1=((tau{chi}(t_chi)'*inv(KI)*tau{chi}(t_chi)))*I/d(chi)/(d(chi)/seno);
Gamma=coeff1*[(1-fx(kk))^2, (1-fx(kk))*fx(kk); fx(kk)*(1-fx(kk)),fx(kk)*fx(kk)];
%Gamma=coeff1*eval(massa);

%aggiungo MEDIA MEDIA
A(Nj(kk):Nj(kk)+1,Nj(kk):Nj(kk)+1)=A(Nj(kk):Nj(kk)+1,Nj(kk):Nj(kk)+1)+0.25*Gamma;
A(Nj(kk):Nj(kk)+1,N+2:N+3)=A(Nj(kk):Nj(kk)+1,N+2:N+3)+0.25*Gamma;
A(extV1:extV2,extV1:extV2)=A(extV1:extV2,extV1:extV2)+0.25*Gamma;
A(extV1:extV2,Nj(kk):Nj(kk)+1)=A(extV1:extV2,Nj(kk):Nj(kk)+1)+0.25*Gamma;

coeff2=(tau{chi}(t_chi)'*inv(KI)*tau{chi}(t_chi))*(d(kk)/(d(chi)));
%Gamma=coeff2*eval(massa);
Gamma=coeff2*[(1-fx(kk))^2, (1-fx(kk))*fx(kk); fx(kk)*(1-fx(kk)),fx(kk)*fx(kk)];

% %aggiungo SALTO SALTO
 segno=[1 -1];
 segno2=[-1 1];
 ind=[Nj(kk) Nj(kk)+1];
 ind2=[extV1 extV2];
 for i=1:2
     for j=1:2
         A(ind(i), ind(j)) = A(ind(i), ind(j)) + segno(i)*segno(j)*0.5*csi0*Gamma(i,j);
         A(ind(i), ind2(j)) = A(ind(i), ind2(j)) + segno(i)*segno2(j)*0.5*csi0*Gamma(i,j);
         A(ind2(i), ind(j)) = A(ind2(i), ind(j)) + segno2(i)*segno(j)*0.5*csi0*Gamma(i,j);
         A(ind2(i), ind2(j)) = A(ind2(i), ind2(j)) + segno2(i)*segno2(j)*0.5*csi0*Gamma(i,j);
     end
 end

 %%aggiungo termine dovuto a fratture non ortogonali


media_media=[1-fx(kk);fx(kk)]*[1-fx_altro(kk), fx_altro(kk)]*(tau{chi}(t_chi)'*inv(KI)*tau{kk}(t_kk)) * I/d(chi)/(d(kk)/seno);
NO{k}=zeros(size(A,1),4);
 segno=[1 -1];
 segno2=[-1 1];
 ind=[Nj(kk) Nj(kk)+1];
 ind2=[extV1 extV2];
 indj=[1 2];
 indj2=[3 4];
%  
 for i=1:2
     for j=1:2
        
          NO{k}(ind(i), indj(j)) = NO{k}(ind(i), indj(j)) + 0.25*media_media(i,j);
         NO{k}(ind(i), indj2(j)) = NO{k}(ind(i), indj2(j)) +0.25*media_media(i,j);
         NO{k}(ind2(i), indj(j)) = NO{k}(ind2(i), indj(j)) + 0.25*media_media(i,j);
         NO{k}(ind2(i), indj2(j)) = NO{k}(ind2(i), indj2(j)) + 0.25*media_media(i,j);
     end
 end


%la C contiene il termine tau valutato nel punto di taglio (per imposizione della pressione)
C(Nj(kk),k)=(1-fx(kk));
C(extV2,k)=fx(kk);
C(Nj(kk)+1,k)=-fx(kk);
C(extV1,k)=-(1-fx(kk));

source(Nj(kk))=source(Nj(kk))*fx(kk);
source(extP)=source(Nj(kk))*(1-fx(kk));

end

F(N+2+2*Nint:end)=source;

%condizioni al contorno

if (t_bc(1)==1)
F(1)=F(1)+p_b(1);
else
F(1)=F(1) - gamma*h(1)^-1*sum(h)*v_bc(1);    
A(1,1)=A(1,1) +gamma*h(1)^-1*sum(h);
end    
if (t_bc(2)==1)
F(N+1)=F(N+1)-p_b(2);
else
F(N+1)=F(N+1) + gamma*h(N)^-1*sum(h)*v_bc(2);    
A(N+1,N+1)=A(N+1,N+1) + gamma*h(N)^-1*sum(h);
end    

