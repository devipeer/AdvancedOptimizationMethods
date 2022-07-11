values = [12, 23, 8, 6, 11, 10, 27, 2];
weights = [2.2, 4, 2.5, 2.1, 1.8, 3, 5, 0.6];
maxWeight = 16;
f = -values;        % max
A = weights;
b = maxWeight;
intcon = 1:8;
lb = zeros(8,1);
ub = ones(8,1);

x = intlinprog(f,intcon,A,b,[],[],lb,ub)