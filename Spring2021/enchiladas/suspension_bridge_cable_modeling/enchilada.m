D = 1;
span = 1;
w = D/span; 
L = sqrt(8);
d = 1;
T = w*L^2/(8*d);

syms y(x)
Dy = diff(y);

ode = diff(y,x,2) == w/T*sqrt(1+diff(y,x)^2);
cond1 = y(0) == 0;
cond2 = Dy(0) == 0;

conds = [cond1 cond2];
ySol(x) = dsolve(ode,conds);

fplot(ySol,[-L/2 L/2])





