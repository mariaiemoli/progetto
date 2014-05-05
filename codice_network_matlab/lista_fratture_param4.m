

%definisco la geometria delle fratture
colori={};

Nhor=5;
Nver=5;
Nf=Nhor+Nver; %numero di fratture


for i=1:Nhor
%F1
colori{i}=[0,0,i/Nhor]*0;
h(i)=1/103;
t{i}=[0:h(i):1];   %nodi
tt{i}=0.5*(t{i}(1:end-1)+t{i}(2:end));  %baricentri
N(i)=length(tt{i});   %numero di elementi
x_map{i}=@(t) [t;  (i-0.5)/(Nhor)+t*0];
Tau{i}=@(t) [1; 0*t];
ds{i}=@(t) (1)^0.5+0*t;
tau{i}=@(t) Tau{i}(t)/ds{i}(t);
%theta(1)=acos()
source{i}=@(t) 0.0*(t<0.7).*(t>0.3)+0*t;
cini{i}=@(t) 0+0*t;
p_bc{i}=[0 0+1*(i==1)]; %condizioni al contorno di pressione
v_bc{i}=0.01*[0 0];
C_bc{i}=[0 0];
x_bc{i}=[0 0];
t_bc{i}=[1 1];
eta(i)=1;  %inverso della permeabilità
mu(i)=100;
d(i)=0.01;%/cos(theta(1)) %spessore
power(i)=1;
xx=x_map{i}(t{i});
plot(xx(1,:),xx(2,:),'-','color',colori{i},'linewidth',2)
hold on
end

for i=Nhor+1:Nver+Nhor
%F1
dth=5*rand(1);
theta=(90+dth)*pi/180;
colori{i}=[(i-Nhor)/Nver,1-(i-Nhor)/Nver,0]*0;
h(i)=1/103;
t{i}=[0:h(i):1];   %nodi
tt{i}=0.5*(t{i}(1:end-1)+t{i}(2:end));  %baricentri
N(i)=length(tt{i});   %numero di elementi
x_map{i}=@(t) [ (i-Nhor-0.5)/(Nver)+t*cos(theta)/sin(theta)-0.5*cos(theta)/sin(theta);t];
Tau{i}=@(t) [ cos(theta)+0*t;sin(theta)];
ds{i}=@(t) (1)^0.5+0*t;
tau{i}=@(t) Tau{i}(t)/ds{i}(t);
%theta(1)=acos()
source{i}=@(t) 0+0*t;
cini{i}=@(t) 0+0*t;
p_bc{i}=[0 1]; %condizioni al contorno di pressione
v_bc{i}=0.01*[0 0];
C_bc{i}=[0 0];
x_bc{i}=[0 0];
t_bc{i}=[0 0 ];
eta(i)=1*(i~=3+Nhor) + 1*(i==3+Nhor);  %inverso della permeabilità
mu(i)=100;
d(i)=0.01;%/cos(theta(1)) %spessore
power(i)=2;
xx=x_map{i}(t{i});
plot(xx(1,:),xx(2,:),'-','color',colori{i},'linewidth',2)
hold on
end


%th=[80 82 85 87.5 90 92.5 95 98 100];
%p=[0.0684 0.0713 0.0761 0.0797 0.0840 0.0886 0.0938 0.1008 0.1058];
%figure
%plot(th,p)
%dens=[3 5 7 9 11 13];
%p=[0.0723 0.0840 0.0901 0.0933 0.0951 0.0961];
%figure
%plot(dens,p)
