%%
clear all
clc

load matrice.mm

%%
M=zeros(max(matrice(:,1)),max(matrice(:,1)));

for i = 1:length(matrice)
    r=matrice(i,1);
    c=matrice(i,2);
    M(r,c)=matrice(i,3);
end

C0=sparse(M);

hold on
spy(C0,'g*');
count=0;

%for i=1:661
%    if M(659,i)~=0
%        count=count+1;
%    end
%end
count