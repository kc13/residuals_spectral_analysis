function [ax] = gridPlot(xmat,ymat,dmat,dtype,suff,figopt)
%gridPlot.m 9/13/22

doYlbl = optCheck(figopt,'doYlbl',[],[]);
doXlbl = optCheck(figopt,'doXlbl',[],[]);
doYRlbl = optCheck(figopt,'doYRlbl',[],[]);
yrmat = optCheck(figopt,'yrmat',[],[]);
hlightCoords = optCheck(figopt,'hlightCoords',[],[]);
prefixMap = optCheck(figopt,'prefixMap',[],[]);

tOut = @(x) set(x,'TickDir','out');

cLimMap = containers.Map;
cLimMap('res') = [0 1];
cLimMap('comp') = [0 1];
cLimMap('diff') = [-0.5 0.5];

if isempty(prefixMap)
    prefixMap = containers.Map;
    prefixMap('res') = 'residuals'; 
    prefixMap('comp') = 'shuffling';
    prefixMap('diff') = [prefixMap('res'),'-',prefixMap('comp')];
end
dstrMap = containers.Map;
dstrMap('hr') = 'hit rates';
dstrMap('fa') = 'false alarm rates';

cLim = cLimMap(suff);

xVals = xmat(1,:);
nX = numel(xVals);
yVals = ymat(:,1);
nY = numel(yVals);

colormap(turbo)

xt = [1:nX]; yt = [1:nY];
im = imagesc(xt,yt,dmat,cLim);
xticks(xt); yticks(yt);
cd = im.CData;
txtstrs = cellstrnum(reshape(cd,[],1));
txfs = 9;
ax = gca;
tOut(ax);
ax.YDir = 'normal';
txtX = reshape(repmat(xt,nY,1),[],1);
txtY = (repmat(yt,1,nX))';
% this is not quite working yet
tx = text(txtX,txtY,txtstrs,'Color',[0.5 0.5 0.5],'FontSize',txfs,...
   'FontWeight','Bold','HorizontalAlignment','Center');
ax = gca;
tOut(ax);
xticklabels(cellstrnum(toCol(xVals)))
yticklabels(cellstrnum(toCol(yVals)))

lblfs = 8;
if doYlbl; ylabel(figopt.ylbl,'FontSize',lblfs); end
if doXlbl; xlabel(figopt.xlbl,'FontSize',lblfs); end

% highlight boxes
if ~isempty(hlightCoords)
    for c = 1:numel(hlightCoords)
        cds = hlightCoords{c}; %x,y pair in param space
        col = find(xVals == cds(1));
        row = find(yVals == cds(2));
        v = [col-0.5 row-0.5; col+0.5 row-0.5; col+0.5 row+0.5; col-0.5 row+0.5];
        f = [1 2 3 4];
        p = patch('Faces',f,'Vertices',v,'FaceColor','none','EdgeColor','k','LineWidth',0.9);  % default lw 0.5
    end
end %if

fig = gcf;

if ~isempty(yrmat)
    yrVals = yrmat(:,1);
    nYR = numel(yrVals);
    yyaxis right
    im = imagesc(xt,yt,dmat,cLim);
    im.AlphaData = zeros(size(dmat));
    ax.YDir = 'normal';
    ytr = [1:nYR];
    yticks(ytr);
    ax = gca;
    tOut(ax);
    yticklabels(cellstrnum(toCol(yrVals)))
    ax.YColor = 'k';
    if doYRlbl; ylabel(figopt.yrlbl,'FontSize',lblfs); end
end %if

tstr = strjoin({prefixMap(suff),dstrMap(dtype)},': ');
title(tstr)

end %fn