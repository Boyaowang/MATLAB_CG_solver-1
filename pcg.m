clc;
clear;
close all;
size =3;

%a = pascal(size);

a=randn(size);
a=a'*a;
B=0.01*eye(size);
C=a+B;

L = randn(size);
L=L'*L;
L = L+B;

b = rand(size,1);
x_real = b'/C;
x0  =rand(size,1);
[x_num,x0_out,i,r] = conjgrad(C,b,x0);
semilogy(r,'o');
hold on;
D = ichol(C);
[x_num,x0_out,i,r] = preconditioned_conjgrad(C,b,x0,L);

semilogy(r);
D= ichol(C);

function [x0, x0_out,i,r_out] = conjgrad(A, b, x0)
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
        if sqrt(rsnew) < 1e-20
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end

end


function a = ichol(a)
	n = size(a,1);
	for k=1:n
		a(k,k) = sqrt(a(k,k));
		for i=(k+1):n
		    if (a(i,k)~=0)
		        a(i,k) = a(i,k)/a(k,k);            
            end
        end
		for j=(k+1):n
		    for i=j:n
		        if (a(i,j)~=0)
		            a(i,j) = a(i,j)-a(i,k)*a(j,k);  
                end
            end
        end
    end

    for i=1:n
        for j=i+1:n
            a(i,j) = 0;
        end
    end            
end

function [x0, x0_out,i,r_out] = preconditioned_conjgrad(A, b, x0,M)
    r = b - A * x0;
    Minv = inv(M);
    %z0 = M\r;
    z0 = Minv * r;
    p = z0;
    rzsold = r' * z0;
    for i = 1:10000 %length(b)
        Ap = A * p;
        alpha = rzsold / (p' * Ap);
        x0 = x0 + alpha * p;
        r = r - alpha * Ap;
        %z = M\r;
        z = Minv * r;
        rzsnew = r' * z;
        x0_out(i,:) = x0;
        r_out(i) = r'*r;
        if sqrt(r'*r) < 1e-20
              break
        end
        p = z + (rzsnew / rzsold) * p;
        rzsold = rzsnew;
    end

end
