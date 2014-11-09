clear all
close all
clc

nel=19192;
%B = ff2matlab('B');
%sigma = ff2matlab('sigma');
%P = ff2matlab('P');

%spy(sigma)
%fopen( 'matriceGrande')
%load ('matriceGrande')
MK =importdata ('matriceGrande');

M=zeros(nel,nel);
%%
for i=1:nel
    for j=1:nel
            r=MK(i,1);
            c=MK(i,2);
            M(r,c)=MK(i,3);
    end
end

MM=sparse(M)

spy(MM)
