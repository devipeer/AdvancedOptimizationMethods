clear all
clc
format short

A = 1;
B = 0.5;
Q = 0.5;
R = 0; % 0 or 0.5
H = 1;
x0 = 1;
N = 6;
K = 25;
eps = 0.25;

u = zeros(K, N);
u(1, :) = [1, 3, 2, 3, 2, 6];
x = zeros(K, N + 1);
x(:, 1) = x0;
p = zeros(1, N + 1);
b = zeros(1, N);
J = zeros(K, 1);
t = 0.01;

for iter = 1 : K
    for i = 2 : size(x, 2)
        x(iter, i) = A * x(iter, i - 1) + B * u(iter, i - 1);
    end
    for i = N : -1 : 1
        p(i) = 2 * x(iter,i) * Q + A * p(i + 1) + R * u(iter,i);
    end
    b_old = b;
    for i = 1 : length(b)
        b(i) = 2 * u(iter, i) * H  + B * p(i + 1) + R * x(iter, i);
    end
    b_norm(i) = sqrt(b*b');
    b_norm_old(i) = sqrt(b_old*b_old');
    if iter == 1
        d = -b;
    else
        d = -b + b_norm(i) ^ 2 / b_norm_old(i) ^ 2 * d;
    end
    J = sum(Q * x(iter, 1 : end - 1).^2 + R*x(iter, end-1)*u(iter, :) + H * u(iter, :).^2);
    J_i(iter) = J;
    if norm(d) < eps
        break;
    else
        syms t;
        x_t = calx(A, B, x0, u(iter, :), -d, t);
        J_t = simplify(calJ(Q,H,x_t,u(iter, :),R, -d,t));
        deriv = diff(J_t);
        T = solve(deriv == 0, t);
        T = double(T);
        u = u + T * d;
    end
end

plot(J_i(1 : iter));
xlabel('Iteration')
ylabel('Performance')
title('Performance index, direct gradient');

function x = calx(A,B,x0,u,b,T)
    N = size(u,2);
    x = sym(zeros(1,N));
    x(1) = sym(x0);
    for i=1:N-1
       x(i+1)=A*x(i)+B*(u(i)-T*b(i)); 
    end
end
function J=calJ(Q, H, x,u,r, b,T)
    N = length(x);
    J = sym(zeros(1,N));
    for i=1:N
       J(i) = (Q*x(i)^2 + r*x(i)*(u(i)-T*b(i)) + H*(u(i)-T*b(i))^2); 
    end
    J = sum(J);
end