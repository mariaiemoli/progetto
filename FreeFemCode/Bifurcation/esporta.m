clear all
close all
clc

%B = ff2matlab('B');
%sigma = ff2matlab('sigma');
%P = ff2matlab('P');

%spy(sigma)
%fopen( 'matriceGrande')
%load ('matriceGrande')
MK =importdata ('matriceGrande');

M=zeros(4829,4829);
%%
for i=1:4829
    for j=1:4829
            r=MK(i,1);
            c=MK(i,2);
            M(r,c)=MK(i,3);
    end
end

MM=sparse(M)

spy(MM)
%%
spy(MK)