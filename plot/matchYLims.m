function [] = matchYLims(goArr,scfrac,minfix,maxfix)
%matchYLims 3/4/21

sz = size(goArr);
if sz(2) > sz(1)
    goArr = goArr';
end

ylArr = cell2mat(arrayfun(@(x) x.YLim,goArr,'UniformOutput',false));
if isempty(minfix)
    ymin = min(ylArr(:,1));
else
    ymin = minfix;
end
if isempty(maxfix)   
    ymax = max(ylArr(:,2));
else
    ymax = maxfix;
end


rg = range([ymin ymax]);
md = median([ymin ymax]);
nrg = rg*(scfrac);
if isempty(maxfix)
    newUB = md + nrg/2; 
else
    newUB = maxfix;
end
if isempty(minfix)
    newLB = md - nrg/2;
else
    newLB = minfix;
end

for g = 1:numel(goArr)
    ax = goArr(g);
    ax.YLim = [newLB newUB];
end

end %fn

