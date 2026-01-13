function nr = exact_norm(A, M, N, b, c, x, y)

MisFunc = isa(M, 'function_handle');
NisFunc = isa(N, 'function_handle');
AisFunc = isa(A, 'function_handle');

if MisFunc; Mx = M(x, 1); else; Mx = M*x; end
if NisFunc; Ny = N(y, 1); else; Ny = N*y; end

if AisFunc
    Ay = A(y, 1);
    ATx = A(x, 2);
else
    Ay = A * y;
    ATx = A' * x;
end

rx = b - (Mx + Ay);
ry = c - (ATx - Ny);

if MisFunc; Mrx = M(rx, 2); else; Mrx = M \ rx; end
if NisFunc; Nry = N(ry, 2); else; Nry = N \ ry; end

nr = sqrt(rx'*Mrx + ry'*Nry);
end