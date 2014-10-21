close all
clear all
clc

data = load('pressure.txt');

x = [];
y = [];
z = [];

index = 1;

for i = data'
    if ( rem(index, 4) ~= 0 )
        x = [x; i(1)+0.5];
        y = [y; i(2)+0.5];
        z = [z; i(3)];
    end
    index = index + 1;
end

tri = [];
nbTriangles = size(x,1) / 3;
for i = 3*[0:nbTriangles-1]
    tri = [tri; [i+1 i+2 i+3]];
end

trisurf(tri,x,y,z,"facecolor","interp")

xlabel('$x$')
ylabel('$y$')
zlabel('$\\primal$')
box('off')
grid('off')
axis([0 2 0 2 -1 1])
pbaspect ([1 1 0.75])
print('-depsc2', '-color', 'pressure.eps')
