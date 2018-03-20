clear;clc;close all
addpath('F:\Matlab Proj\GeneralUtility');
% r = 1;
% R = 10+r;
% position = importfile('position.dat');
% n = length(position);
% fai = n*r*r/R/R;
% for i = 1:n
%     PlotSolidCircle(position(i,1), position(i,2), r )
% end
% axis equal
% 
% 
% Fai = 0.5;
% a = R*sqrt(pi/Fai);
% axis([-a,a, -a,a]);
% set(gca,'xtick',[],'xticklabel',[])
% set(gca,'ytick',[],'yticklabel',[])
% set(gca, 'box','on')
% set(gca,'position',[0 0 1 1])
% set(gcf,'position',[0 0 788 788])
% saveas(gca, 'Gen.pgm', 'pgm')
A = imread('Gen.pgm');
A(end,:) = 255;
A(:,1) = 255;
A(:,end) = 255;
A(1,:) = 255;
imshow(A)
A=imresize(A,1.5);
figure
imshow(A)

A(A~=255)  = 1;
A(A==255)  = 0;

[Nx,Ny] = size(A);
A = [A;zeros(1000,Ny)];
[Nx, Ny] = size(A);
A_ = reshape(A, 1, Nx*Ny);
dlmwrite('PPGeo.dat', A_, 'delimiter',' ');




