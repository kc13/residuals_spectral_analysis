function [xCrop,yCrop] = cropToOneDimRange(x,y,xLB,xUB)
%cropToOneDimRange.m

    assert(isvector(x) && isvector(y),'only accepts vector inputs currently')

    LBleftIX = find(x <= xLB,1,'last');
    LBrightIX = find(x >= xLB,1,'first');
    
    % LB
    if isempty(find(x == xLB))
        flankX = [x(LBleftIX) x(LBrightIX)];
        flankY = [y(LBleftIX) y(LBrightIX)];
        yLB = interp1(flankX,flankY,xLB);
        ix = LBleftIX+1;
        xExpand = insertValue(x,xLB,ix);
        yExpand = insertValue(y,yLB,ix);
    else
        xExpand = x;
        yExpand = y;
    end
    
    UBleftIX = find(xExpand <= xUB,1,'last');
    UBrightIX = find(xExpand >= xUB,1,'first');  
    % UB
    if isempty(find(xExpand == xUB))
        flankX = [xExpand(UBleftIX) xExpand(UBrightIX)];
        flankY = [yExpand(UBleftIX) yExpand(UBrightIX)];
        yUB = interp1(flankX,flankY,xUB);
        ix = UBleftIX+1;
        xExpand2 = insertValue(xExpand,xUB,ix);
        yExpand2 = insertValue(yExpand,yUB,ix);
    else
        xExpand2 = xExpand;
        yExpand2 = yExpand;
    end
    
    ix = in_range_incl(xExpand2,[xLB,xUB]);
    xCrop = xExpand2(ix);
    yCrop = yExpand2(ix);
    
    assert(min(xCrop) == xLB && max(xCrop) == xUB && issorted(x),...
        'cropped x vector not as expected')
    try
    assert(issorted(yCrop) && nnz(~ismember(yCrop,y)) <= 2,...
        'cropped y vector not as expected')  
    catch
		disp('unexpected behavior in cropToOneDimRange: debug?')
        keyboard
    end    
    
end
