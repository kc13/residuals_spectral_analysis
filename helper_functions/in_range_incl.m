function [y] = in_range_incl(x,rnge)

    y = x >= rnge(1) & x <= rnge(2);

end

