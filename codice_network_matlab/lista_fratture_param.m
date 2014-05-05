
Nf=4; %numero di fratture

%definisco la geometria delle fratture
colori={};


%F1
colori{1}=[0,0,1];
h(1)=1/103;
t{1}=[0:h(1):1];   %nodi
tt{1}=0.5*(t{1}(1:end-1)+t{1}(2:end));  %baricentri
N(1)=length(tt{1});   %numero di elementi
x_map{1}=@(t) [t-0.5;  2*t-1];
Tau{1}=@(t) [1; 2];
ds{1}=@(t) (5)^0.5+0*t;
tau{1}=@(t) Tau{1}(t)/ds{1}(t);
%theta(1)=acos()
source{1}=@(t) 0+0*t;
cini{1}=@(t) 0+0*t;
p_bc{1}=[0 1]; %condizioni al contorno di pressione
v_bc{1}=0.01*[-1 0];
C_bc{1}=[0 0];
x_bc{1}=[0 0];
t_bc{1}=[1 1];
eta(1)=1;  %inverso della permeabilità
mu(1)=100;
d(1)=0.005;%/cos(theta(1)) %spessore
power(1)=1;
xx=x_map{1}(t{1});
plot(xx(1,:),xx(2,:),'-o','color',colori{1},'linewidth',2)
hold on


%F2
colori{2}=[1,0,0];
h(2)=1/105;
t{2}=[0:h(2):1];   %nodi
tt{2}=0.5*(t{2}(1:end-1)+t{2}(2:end));  %baricentri
N(2)=length(tt{2});   %numero di elementi
x_map{2}=@(t) [t-0.5;  -2*t+1.1];
Tau{2}=@(t) [1; -2];
ds{2}=@(t) (5)^0.5+0*t;
tau{2}=@(t) Tau{2}(t)/ds{2}(t);
%theta(1)=acos()
source{2}=@(t) 0+0*t;
cini{2}=@(t) 0+0*t;
p_bc{2}=[0 1]; %condizioni al contorno di pressione
v_bc{2}=[0 0];
C_bc{2}=[1 0];
x_bc{2}=[0 0];

t_bc{2}=[1 1];
eta(2)=100;  %inverso della permeabilità
mu(2)=100;
d(2)=0.005;%
power(2)=2;
xx=x_map{2}(t{2});
plot(xx(1,:),xx(2,:),'-o','color',colori{2},'linewidth',2)
hold on

% %F3
colori{3}=[0,0,0];
h(3)=1/105;
t{3}=[0:h(3):1];   %nodi
tt{3}=0.5*(t{3}(1:end-1)+t{3}(2:end));  %baricentri
N(3)=length(tt{3});   %numero di elementi
x_map{3}=@(t) [(t-0.5);  sin(3*t)];
Tau{3}=@(t) [1; 3*cos(3*t)];
ds{3}=@(t) (1+9*cos(3*t).^2).^0.5;
tau{3}=@(t) Tau{3}(t)/ds{3}(t);
%theta(1)=acos()
source{3}=@(t) 0+0*t;
cini{3}=@(t) 0+0*t;
p_bc{3}=[0 1]; %condizioni al contorno di pressione
v_bc{3}=[0 0];
C_bc{3}=[0 0];
x_bc{3}=[0 0];

t_bc{3}=[1 1];
eta(3)=1;  %inverso della permeabilità
mu(3)=100;
d(3)=0.005;%
power(3)=1;
xx=x_map{3}(t{3});
plot(xx(1,:),xx(2,:),'-o','color',colori{3},'linewidth',2)
hold on
% 
%F3
colori{4}=[0,0.6,0];
h(4)=1/105;
t{4}=[0:h(4):1];   %nodi
tt{4}=0.5*(t{4}(1:end-1)+t{4}(2:end));  %baricentri
N(4)=length(tt{4});   %numero di elementi
x_map{4}=@(t) [(t-0.5);  1.2*t-0.85];
Tau{4}=@(t) [1; 1.2+t*0];
ds{4}=@(t) (1+1.2^2+0*t).^0.5;
tau{4}=@(t) Tau{4}(t)/ds{4}(t);
%theta(1)=acos()
source{4}=@(t) 0+0*t;
cini{4}=@(t) 0+0*t;
p_bc{4}=[0 1]; %condizioni al contorno di pressione
v_bc{4}=[0 0];
C_bc{4}=[0 1];
x_bc{4}=[0 0];

t_bc{4}=[1 1];
eta(4)=1;  %inverso della permeabilità
mu(4)=100;
d(4)=0.005;%
power(4)=1;
xx=x_map{4}(t{4});
plot(xx(1,:),xx(2,:),'-o','color',colori{4},'linewidth',2)
% hold on
% 
% 
% %F5
% colori{5}=[1,0.0,0.8];
% h(5)=1/67;
% t{5}=[0:h(5):1];   %nodi
% tt{5}=0.5*(t{5}(1:end-1)+t{5}(2:end));  %baricentri
% N(5)=length(tt{5});   %numero di elementi
% x_map{5}=@(t) [-0.15+0*t; 2.5*t-1.5];
% Tau{5}=@(t) [0; 2.5+0*t];
% ds{5}=@(t) 2.5+0*t;
% tau{5}=@(t) Tau{5}(t)/ds{5}(t);
% %theta(1)=acos()
% source{5}=@(t) 0+0*t;
% cini{5}=@(t) 0+0*t;
% p_bc{5}=[0 0]; %condizioni al contorno di pressione
% v_bc{5}=-0.005*[1 1];
% C_bc{5}=[0 0];
% x_bc{5}=[0 0];
% 
% t_bc{5}=[1 1];
% eta(5)=1;  %inverso della permeabilità
% mu(5)=100;
% d(5)=0.005;%
% power(5)=1;
% xx=x_map{5}(t{5});
% plot(xx(1,:),xx(2,:),'-o','color',colori{5},'linewidth',2)
% hold on
% 
% 
% %F6
% colori{6}=[0,0.8,0.8];
% h(6)=1/67;
% t{6}=[0:h(6):1];   %nodi
% tt{6}=0.5*(t{6}(1:end-1)+t{6}(2:end));  %baricentri
% N(6)=length(tt{6});   %numero di elementi
% x_map{6}=@(t) [t-0.5; -t+0.7];
% Tau{6}=@(t) [1; -1+0*t];
% ds{6}=@(t) (2)^0.5+0*t;
% tau{6}=@(t) Tau{6}(t)/ds{6}(t);
% %theta(1)=acos()
% source{6}=@(t) 0+0*t;
% cini{6}=@(t) 0+0*t;
% p_bc{6}=[0 0 ]; %condizioni al contorno di pressione
% v_bc{6}=[-1 -1];
% C_bc{6}=[0 0];
% x_bc{6}=[0 0];
% t_bc{6}=[1 1];
% 
% eta(6)=1;  %inverso della permeabilità
% mu(6)=100;
% d(6)=0.005;%
% power(6)=1;
% xx=x_map{6}(t{6});
% plot(xx(1,:),xx(2,:),'-o','color',colori{6},'linewidth',2)
% hold on
% 
