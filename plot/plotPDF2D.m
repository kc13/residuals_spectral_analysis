function [ax] = plotPDF2D(x,y,data,figprefs)
%plotPDF2D.m
    
    cmap = colormap(turbo);
    tickdir = 'out';
    dobox = 'off';
    
    xL = optCheck(figprefs,'xlims',[],[]);
    yL = optCheck(figprefs,'ylims',[],[]);
    xLbl = optCheck(figprefs,'xlbl','ISI(n) (ms)',[]);
    yLbl = optCheck(figprefs,'ylbl','ISI(n+1) (ms)',[]);
    clims = optCheck(figprefs,'clims',[],[]);

    im = imagesc(x,y,data);
    colorbar 
    ax = gca;
    ax.TickDir = tickdir;
    ax.Box = dobox;
    ax.YDir = 'normal';

    im.AlphaData = ones(size(im.CData));
    im.AlphaData(im.CData == 0) = 0;

    if ~isempty(xL)
        ax.XLim = xL;
    end
    if ~isempty(yL)
        ax.YLim = yL;
    end
    if ~isempty(clims)
        ax.CLim = clims;
    end
    xlabel(xLbl)
    ylabel(yLbl)

    title('joint PDF: ISI(n), ISI(n+1)')


end