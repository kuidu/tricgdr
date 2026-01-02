function y = Mopertor(M, dM, x, k)
if k == 1
    y = M * x;
else
    y = dM \ x;
end
end
