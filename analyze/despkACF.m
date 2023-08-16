function [acf2] = despkACF(acf)
%despkACF.m 4/25/22
acf2 = acf;
acf2(ceil(length(acf)/2)) = nan;
end

