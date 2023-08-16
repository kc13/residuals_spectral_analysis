function [tl] = tlcompact(nR,nC,titlestr,tftsz)
%tlcompact.m
%title string and title font size are optional

tl = tiledlayout(nR,nC,'Padding','compact','TileSpacing','compact');

n = nargin;
if n >= 3 
    tl.Title.String = titlestr;
    if n > 3
        tl.Title.FontSize = tftsz;
    end
end


end

