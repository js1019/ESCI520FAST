function power = specCalc(s)
for n = 1:size(s,2)
    for m = 1:size(s,1)
        power(m,n) = real(s(m,n))^2+imag(s(m,n))^2;
    end
end
end
