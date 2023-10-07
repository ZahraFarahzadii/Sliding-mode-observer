%%  ALLAH

function dX = odeEquations(t,X)
global A B F C Gn Gl
global normey
xhat = X(1:4,1);
x = X(5:8,1);
f = 0.6*sin(5*t);
u = -lqr(A,B,diag([1,1,1,1]),1)*x;
ey = C*xhat - C*x;
if t == 0
    normey = 0;
else
    normey = norm([normey;ey]);
end
v = -ey/(normey + eps);
dxhat = A*xhat + B*u - Gl*ey + Gn*v;
dx = A*x + B*u + F*f;
dX = [dxhat;dx];
end