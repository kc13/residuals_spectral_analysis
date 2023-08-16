function [HistSpl] = getSplines(method,lag,c_pt,s)
% calls functions from
%https://github.com/MehradSm/Modified-Spline-Regression
% ref:
%Sarmashghi, M., Jadhav, S. P., & Eden, U. (2021). 
%Efficient spline regression for neural spiking data. 
%Plos one, 16(10), e0258321.

switch method
    case 'modcard'
        if numel(c_pt) == 2
            [HistSpl] = OneModifiedCardinalSpline(lag,c_pt);
        else
            [HistSpl] = ModifiedCardinalSpline(lag,c_pt,s);
        end
    case 'card'
        if numel(c_pt) < 4
            error('cardinal splines require >= 4 control points')
        else
            [HistSpl] = CardinalSpline(lag,c_pt,s);
        end
    otherwise
        error('unrecognized spline method')
end %sw
            

end
