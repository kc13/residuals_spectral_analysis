function [HistSpl] = OneModifiedCardinalSpline(lag,c_pt)
%6/19/22
% adapting method from Sarmashghi et al. paper for one MC spline

assert(numel(c_pt) == 2,'this function only handles two control points')

HistSpl = zeros(lag,length(c_pt));
%for each 1 ms timepoint, calculate the corresponding row of the glm input matrix
for i = 1:lag
    nearest_c_pt_index = max(find(c_pt<i)); 
    nearest_c_pt_time = c_pt(nearest_c_pt_index);
    next_c_pt_time = c_pt(nearest_c_pt_index+1);
    % Compute the fractional distance between timepoint i and the nearest knot
    u = (i-nearest_c_pt_time)./(next_c_pt_time-nearest_c_pt_time);
    p=[u^3 u^2 1]*[2 2; -3 3; 1 0];
    HistSpl(i,:) = p;
end %i



end

% original spline refs:
%https://github.com/MehradSm/Modified-Spline-Regression
% ref:
%Sarmashghi, M., Jadhav, S. P., & Eden, U. (2021). 
%Efficient spline regression for neural spiking data. 
%Plos one, 16(10), e0258321.