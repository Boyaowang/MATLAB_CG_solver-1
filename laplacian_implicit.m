clear;
clc;
close all;

dx = 1;
dy =1;
dt =1;
t0 = dt;
T=5;
A = dx*dy/dt;
B= dy/dx;

M = zeros(9,9);

M(1,1)= B+1/B;
M(3,3) = M(1,1);
M(7,7) =M(1,1);
M(9,9) = M(1,1);
M(2,2) = B+2/B;
M(8,8) = M(2,2);
M(4,4) =  2*B+1/B;
M(6,6) =M(4,4);

M(5,5) = 2*B +2/B;

for i = 1:9
    M(i,i) = M(i,i)+ dt;
end

for i = 1:6
    M(i,i+3) = -B;
    M(i+3, i)=-B;
end

for i=1:9
    if(i ==1 || i ==2 || i ==4 ||i ==5 ||i ==7 ||i ==8)
        M(i,i+1) = -1/B;
    end
    
    if(i ==2 || i ==3 || i ==5 ||i ==6 ||i ==8 ||i ==9)
        M(i,i-1) = -1/B;
    end
end

S = [0 0 0 0 1 0 0 0 0]';

b0 = S/dt;
S_save = zeros(9,(T/dt));

for i= 1:T/dt
    [S,r_out] = conjgrad(M, b0, S);
    S_save(:,i) = S;
    b0 = S/dt;
end

TempF = [1 2 4 1 5 1 4 1 4]';

for i= 1:T/dt
    S = dt * M * S_save(:,T+1-i) .* TempF;
end


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


        