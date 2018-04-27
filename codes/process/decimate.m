function [dt, series, tNew] = decimate(dec, series, tc)
dt = 1/dec;
dtLength = dt/(tc(2)-tc(1));
count = 1;
for n = 1:length(series)
    if rem(n,int32(dtLength)) == 0
        temp(count) = series(n);
        tNew(count) = tc(n);
        count = count + 1;
    end 
end
clear series
series = temp;
end