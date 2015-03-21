clc;
clear; 
randn('seed',1)

n = 2^11;
m = 2^11;
block_size = 2^7;
A = randn(m,n);

machine_pre =n*eps(1)

%BlockQR
disp('BlockQR')
tic
[Q R] = blockQR(A, block_size);
time = toc

Gflop = (2*m*n*n)/(time*1e9)
max(max(Q'*Q -eye(n))) 
max(max(A-Q*R)) 

%BlockReflect
disp('Block Reflect')
tic
[Q R] = blockReflect(A);
time = toc

Gflop = (2*m*n*n)/(time*1e9)
max(max(Q'*Q -eye(n))) 
max(max(A-Q*R)) 

%MATLAB QR
disp('MATLAB QR')
tic
[Q R] = qr(A);
time = toc
Gflop = (2*n*n*n)/(time*1e9)
max(max(Q'*Q -eye(n))) 
max(max(A-Q*R))