clear;
clc;
close all;


%% Start descritization %%
L = 0.09;
dx = 0.01;
dt = 0.0001;
T = 0.001;
D = 0.001;

alpha = dx/dt;
beta = D/dx;

N = L/dx;
M = zeros(N,N);

for i = 2:N-1
    M(i,i) = alpha+2*beta;
    M(i,i-1) = -beta; M(i,i+1) = -beta;
end

M(1,1) = alpha+beta; M(1,2) = -beta;
M(N,N) = alpha+beta; M(N,N-1) = -beta;

%% Finish descritization  and do the calculation %%

S = [1 0 0 0 1 0 0 0 0]';
b0 = alpha.*S;
S_save = zeros(N,(T/dt));
S0 = S;

%% diffusion of S %%
for i=1:T/dt
    [S,r_out] = conjgrad(M, b0, S);
    S_save(:,i) = S;
    b0 = S*alpha;
end


TempF = [1 0 6 2 6 2 6 0 6]';
STp = S' * TempF
 

%% diffusion of T %%
 bt = alpha .*  TempF;
 
 for i=1:T/dt
    [TempF,r_out] = conjgrad(M, bt, TempF);
    bt = TempF*alpha;
 end
 
%% processing T %%

ST_par = S0' * TempF
ST_tempDiffusion = S0 .* TempF
plot(ST_tempDiffusion, '-o')
% ST = (ST./sum(ST));
% ST = ST .* STp

%  
%  figure(2)
% plot(ST_par)

sum(ST_tempDiffusion)



%% The CG method for solving the linear system %%

function [x0,r_out] = conjgrad(A, b, x0)
    r = b - A * x0;
    p = r;
    rsold = r' * r;
    for i = 1:10000 %length(b)
        Ap = A * p;
        alpha = 1 * (rsold / (p' * Ap));
        x0 = x0 + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        x0_out(i,:) = x0;
        r_out(i) = rsnew;
        if sqrt(rsnew) < 1e-5
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end

end