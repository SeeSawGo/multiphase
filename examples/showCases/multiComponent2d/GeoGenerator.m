clear;clc
Nx = 400;
Ny = 100;
width1 = 50;
width2 = 20;
A = zeros(Ny, Nx, 'int8');

A(width1:Ny, 1:100) = 1;
A(width1:Ny, 100+width2:Nx) = 1;

A(1, :) = 1; %bottom
%A(Ny,:) = 1;
% A(:, 1) = 1;
% A(:, Nx) = 1;

imagesc(A)
A_ = reshape(A, 1, Nx*Ny);
dlmwrite('Geo.dat', A_, 'delimiter',' ');