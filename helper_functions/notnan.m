function [cleanvec] = notnan(vec)

    cleanvec = vec(~isnan(vec));

end %fn