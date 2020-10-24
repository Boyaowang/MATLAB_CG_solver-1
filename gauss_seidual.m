clc;
clear;
close all;
size =100;

%a = pascal(size);

a=randn(size);
a=a'*a;
B=0.01*eye(size);
C=a+B;

b = rand(size,1);
x_real = b'/C;
x0  =rand(size,1);
[x_num,x0_out,i,r] = conjgrad(C,b,x0);
semilogy(r);

function [x0, x0_out,i,r_out] = conjgrad(A, b, x0)
    r = b - A * x0;
    %p = r;
    p= zeros(length(b),1);
    rsold = r' * r;
    for i = 1:length(b)
        p(i,1)=1;
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x0 = x0 + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        x0_out(i,:) = x0;
        r_out(i) = rsnew;
        if sqrt(rsnew) < 1e-20
              break
        end
        %p = r + (rsnew / rsold) * p;
        p= zeros(length(b),1);
        rsold = rsnew;
    end

end