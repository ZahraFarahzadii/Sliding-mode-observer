%%  ALLAH

function ys = lagrangePolynomial(t,y,ts)
ys = 0;
for i = 1:length(t)
    Y1 = 1;
    for j = 1:length(t)
        if i ~= j
            Y1 = Y1*(ts - t(j))/(t(i) - t(j));
        end
    end
    Y1 = y(i)*Y1;
    ys = ys + Y1;
end
end