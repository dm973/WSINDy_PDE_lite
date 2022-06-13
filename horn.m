function y=horn(c,x)
    l = length(c);
    y = c(1);
    for j=1:l-1
        y = y.*x+c(j+1);
    end
end