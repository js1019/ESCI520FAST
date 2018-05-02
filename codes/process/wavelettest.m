% apply wavelet transform
clear all; clc;
load('../../data/sSample.mat');
x = log(sSample);

[nx,ny] = size(x);
% 2D Haar wavelet transform; 
% it is introduced in matlab r2016b
%[a,h,v,d] = haart2(sSample); 

% 2D inverse Haar wavelet transform
%xrec = ihaart2(a,h,v,d);

% wavelet decomposition
% wavelet level
N = 2;
[c0,s]=wavedec2(x,N,'haar');
%figure; imagesc(reshape(c0,nx,ny));
% reconstruction
xrec0 = waverec2(c0,s,'haar');

% reconstruction error
norm(xrec0-x)/norm(x)


c1 = c0;
[cn,ord] = sort(abs(c0)); 
num = 800; 
% play with threshold
thred = abs(c0(ord(nx*ny-num)));
c1(ord(1:end-num)) = 0; 
%figure; imagesc(reshape(c1,nx,ny));
xrec1 = waverec2(c1,s,'haar');

% reconstruction error
norm(xrec1-x)/norm(x)

c2 = c1; 
c2(c0>thred) = 1; 
c2(c0<-thred) = -1;
% number of nonzeros
sum(abs(c2)==1)
%**********************************
% c2: modified wavelet coefficients
%**********************************
figure; imagesc(reshape(c2,nx,ny));

cnew0 = reshape(c2,nx,ny);
cnew1 = zeros(nx*2,ny); 
ct = zeros(nx,ny); cs = zeros(nx,ny);

ct(cnew0<0.5) = 0; 
ct(cnew0>0.5) = 1; 
cnew1(1:2:nx*2-1,:) = ct;
cs(cnew0>-0.5) = 0;
cs(cnew0==-1) = 1; 
cnew1(2:2:nx*2,:) = cs;
figure; imagesc(cnew1)

