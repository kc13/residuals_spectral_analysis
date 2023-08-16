function [ax] = plotMdlDD(rpInfo,HLix,priColors,figprefs)
%plotMdlDD.m

dd = cellfun(@(x) x.chi2Stat(2),rpInfo.dt);
optLag = rpInfo.optLag;
isiLags = rpInfo.isiLags;
ddLags = isiLags(1:end-1);

tickdir = 'out';
dobox = 'off';
dogrid = false;

doarrow = optCheck(figprefs,'doarrow',false,[true false]);
doleg = optCheck(figprefs,'doleg',true,[false true]);

plot(ddLags,dd,'Color',priColors('bars'))
hold all;

HLclrs = priColors('HL');
nHL = numel(HLix);
goArr = gobjects(nHL,1);
lw = 1; %0.5 default
mksm = 6; % also default
mklg = 18;
mkopt = 10;

nL = numel(ddLags);
for L = 1:nL
    if ~ismember(L,HLix)
        mksz = mksm;
        if isiLags(L) <= optLag+1
            mkclr = 'k';
        else
            mkclr = priColors('bars');
        end
        plot(L,dd(L),'.','MarkerSize',mksz,'Color',mkclr) 
    else
        ix = find(HLix == L);
        if isiLags(L) == optLag
            mksz = mkopt;
            mksh = 'p';
        else
            mksz = mklg;
            mksh = '.';
        end
        mkclr = HLclrs(ix,:);
        goArr(ix) = plot(L,dd(L),mksh,'MarkerSize',mksz,'Color',mkclr,'MarkerFaceColor',mkclr);
    end
end %L

hold off;
if dogrid; grid on; end

xlabel('exponential fit offsets (ms)')
ylbl = optCheck(figprefs,'ylbl','D(constant)-D(constant+lag)',[]);
ylabel(ylbl)

if doleg
    legtext = regexprep(cellstrnum(toCol(HLix)),'(.+)','start lag = $1 ms');  
    legend(goArr,legtext,'Box','off','location','Southeast')
end

if isfield(figprefs,'xlims')
   xlim(figprefs.xlims) 
end

ax = gca;
yL = ax.YLim;

if doarrow
    ta = text(rpInfo.RPend,dd(rpInfo.RPend)-0.05*diff(yL),'\uparrow','FontSize',13,...
         'HorizontalAlignment','center','VerticalAlignment','top'); % default units data
    tt = text(rpInfo.RPend,dd(rpInfo.RPend)-0.25*diff(yL),'RP end','FontSize',ax.FontSize,...
        'HorizontalAlignment','center','VerticalAlignment','top');
end %if

doTri = optCheck(figprefs,'doTri',[false],[true false]);
if doTri
    hold on;
    t = plot(rpInfo.RPend,dd(rpInfo.RPend),'^b');
    ax.Parent.Children(1).String{2} = 'RP end';
end

ax.TickDir = tickdir;
ax.Box = dobox;

tstr = 'exponential fit: deviance test';
title(tstr)

end %fn

