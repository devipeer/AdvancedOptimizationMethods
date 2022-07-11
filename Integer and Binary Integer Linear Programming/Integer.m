ObjF = [2.7572 5.5897];
Cons = [-1.9785 4.6353 5.5897 2.2886 2.7572 6.6232];
f = -1 * ObjF;
A = [-1 0; 0 -1; Cons(1,4)/Cons(1,1) 1; Cons(1,6)/Cons(1,2) 1; Cons(1,5)/Cons(1,3) 1];
b = [0; 0; Cons(1,4); Cons(1,6); Cons(1,5)];

options = optimoptions('linprog','Display','none');
[x, fval] = linprog(f, A, b, [], [], [], [], options);
% x2 <= 2 and x2 >= 3
tmp1 = [0 1];
A1 = [A; tmp1];
b1 = [b; floor(x(2))];
A2 = [A; -1*tmp1];
b2 = [b; -ceil(x(2))];
[roz1, fval1] = linprog(f, A1, b1, [], [], [], [], options);
[roz2, fval2] = linprog(f, A2, b2, [], [], [], [], options);

% x1 <= 1 and x1 >= 1
tmp1 = [1 0];
A1 = [A1; tmp1];
b1 = [b1; floor(roz1(1))];
A2 = [A2; -1*tmp1];
b2 = [b2; -ceil(roz1(1))];
[roz1, fval1] = linprog(f, A1, b1, [], [], [], [], options);
[roz2, fval2] = linprog(f, A2, b2, [], [], [], [], options);