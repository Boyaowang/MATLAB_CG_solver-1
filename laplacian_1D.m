clear;
clc;
close all;

L = 0.09;
dx = 0.01;
dt = 0.0001;
T = 0.001;
D = 0.01;

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


S = [1 0 0 0 3 0 3 0 1]';
b0 = alpha.*S;
S_save = zeros(N,(T/dt));

for i=1:T/dt
    [S,r_out] = conjgrad(M, b0, S);
        S_save(:,i) = S;
    b0 = S*alpha;
end


TempF = [4 1 3 7 4 7 3 1 4]';

ST = S.*TempF


sum(ST)

 for i= 1:(T/dt)+10
     ST =(M * ST)./alpha;
     
     for ii = 1:50
     for i = 1:N
        if(ST(i) <0)
              if(i==1)
                  ST(i+1) = ST(i+1)+ST(i)/2;
                  ST(i+2) = ST(i+2)+ST(i)/2;
                  ST(i) = 0;
              elseif(i==N)
                  ST(i-1) = ST(i-1)+ST(i)/2;
                  ST(i-2) = ST(i-2)+ST(i)/2;
                  ST(i) = 0;
              else
              ST(i-1) = ST(i-1)+ST(i)/2;
              ST(i+1) = ST(i+1)+ST(i)/2;
              ST(i) = 0;
              end
        end
     end
     end
     plot(ST)
     hold on
 end
 
 figure(2)
 plot(ST)
 
%  bst = alpha.*ST;
 
%      plot(ST);
%      hold on;

%  for i=1:20
%     [ST,r_out] = conjgrad(M, bst, ST);
%     bst = ST*alpha;
%     S_save(:,i) = S;
%     plot(ST);
%     hold on;
%  end
%  
%   for i= 1:10
%      ST =(M * ST)./alpha;
%  end
 
sum(ST)

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