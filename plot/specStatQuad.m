function [fig] = specStatQuad(xcData,sData,rData,quadopt)


    oti = optCheck(quadopt,'oti','acf + original, corrected PSDs',[]);
    XL = optCheck(quadopt,'XL',[0 500],[]);
    XTdt = optCheck(quadopt,'XTdt',50,[]);
    axftsz = optCheck(quadopt,'axftsz',[],[]);
    figopt.axftsz = axftsz;
    lw = optCheck(quadopt,'lw',[],[]);
    figopt.lw = lw;

    % argument checks
    flbl = optCheck(quadopt,'flbl','frequency (Hz)',[]);
    parent = optCheck(quadopt,'parent',[],[]);
    tilenum = optCheck(quadopt,'tilenum',[],[]);
    plotShuf = optCheck(quadopt,'plotShuf',false,[true false]);
    legloc = optCheck(quadopt,'legloc','Southeast',[]);
    otifs = optCheck(quadopt,'otifs',10.5,[]);

    % check if delta functions or splines
    resfld = optCheck(quadopt,'resfld','res',{'res','sRes'});

    if isempty(parent)
	    TL = tlcompact(2,2,oti,otifs); 
    else
        TL = tiledlayout(parent,2,2,'Padding','compact','TileSpacing','compact'); 
        TL.Layout.Tile = tilenum;
        TL.Title.String = oti;
        TL.Title.FontSize = otifs;
    end
	
    nT = prod(TL.GridSize);
	dataArr = {xcData,sData,sData,rData};
	f1Arr = {[],[],'comp',resfld};
	f2Arr = {[],'S1','S1Comp','S1'};
    tstrArr = {[],'uncorrected','shuffling','residuals'};
    
	for t = 1:nT
		ct = nexttile(TL,t);
        data = dataArr{t};
        f1 = f1Arr{t};
        f2 = f2Arr{t};
        tstr = tstrArr{t};
        figopt = struct;
		if t == 1
            figopt.ax = ct;
            RPend = rData.(resfld).rpInfoS.RPend; 
            if isstruct(rData.(resfld).rpInfoS.RPend)
               figopt.tistr = sprintf('ACF: mean FR = %0.2f Hz, burst, nonburst RP = %u, %u ms',...
                  data.src_FR,RPend.b,RPend.nb);
            else
                figopt.tistr = sprintf('ACF: mean FR = %0.2f Hz, estimated RP = %u ms',...
                    data.src_FR,rData.(resfld).rpInfoS.RPend);
            end
            acorrFig(data,figopt);
        else
            if t == 2
                figopt.doleg = true;
                figopt.legloc = legloc;
                figopt.plotShuf = plotShuf;                  
            else
                figopt.doleg = false;
                figopt.plotShuf = false;
            end
            figopt.ax = ct;
            figopt.xlim = XL; % why is this here also
            figopt.XTdt = XTdt;           
            figopt.tstr = tstr;
            figopt.XL = XL;
            figopt.flbl = flbl;
		    specStatFig(data,f1,f2,figopt);
		end % if

    ax = gca; 
    fig = gcf; 
    
end %fn