clc;
clear;
close all;
size =10;
a = pascal(size);
b = rand(size,1);
x_real = b'/a;
x0  =rand(size,1);
[x_num,x0_out,i,r] = conjgrad(a,b,x0);
semilogy(r);

function [x0, x0_out,i,r_out] = conjgrad(A, b, x0)
    r = b - A * x0;
    p = r;
    rsold = r' * r;
    for i = 1:1000 %length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x0 = x0 + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        x0_out(i,:) = x0;
        r_out(i) = rsnew;
        if sqrt(rsnew) < 1e-10
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end

end
