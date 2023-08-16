function [ax] = specStatFig(specStruct,fld1,fld2,figopt)
%specStatFig.m 

	clrsMap = containers.Map();
	clrsMap('sig') = [0 1 0];
	clrsMap('orig') = [0 0 0];
	% comp = shuffling "compensated" (orig/mean shuf)
	clrsMap('comp') = [0.3020 0.0784 0.5490]; 
	clrsMap('thr') = ones(1,3)*0.35;
	clrsMap('shuf') = [1 0 0];
	clrsMap('res') = [0.0039 0.1255 0.2353];
	clrsMap('sRes') = clrsMap('res'); % sRes is spline output
    Lmap = containers.Map();
	Lmap('orig') = 'original';
	Lmap('comp') = 'original/shuffled';
	Lmap('thr') = '95% CI (corrected, one-sided)';
    Lmap('sig') = 'p < 0.05 & \in (0 100]';
    Lmap('shuf') = 'mean(shuffled)';
    Lmap('res') = 'residuals';
    Lmap('resLag') = Lmap('res');
    Lmap('sRes') = 'residuals (splines)';
    % default LineWidth: 0.5000
    Wmap = containers.Map();
    Wmap('thr') = 0.9;

    XTdt = optCheck(figopt,'XTdt',50,[]);
    XT = optCheck(figopt,'xticks',[],[]); % will override XTdt if set
    if ~isempty(XT)
        XTLmod = median(diff(XT));
    else
        XTLmod = XTdt;
    end
    flbl = optCheck(figopt,'flbl','frequency (Hz)',[]);	
	
    axftsz = optCheck(figopt,'axftsz',[],[]);
    lw = optCheck(figopt,'lw',[],[]);
	
	tickdir = 'out';
	
    dogrid = false;
    dobox = 'off';	
	
    mksz = 20;
    
    if ~exist('figopt','var')
        figopt = struct();
    end	

    arrowfreq = optCheck(figopt,'arrowfreq',[],[]);

	doleg = optCheck(figopt,'doleg',true,[true false]);
    if doleg
        legloc = optCheck(figopt,'legloc','Northeast',{'Northeast','Northwest','Southeast','Southwest','Best'});
    end
	plotShuf = optCheck(figopt,'plotShuf',false,[true false]);
    
    inputax = optCheck(figopt,'ax',[],[]);
    XL = optCheck(figopt,'xlim',[0 500],[]);
    YL = optCheck(figopt,'ylim',[],[]);
    XT = optCheck(figopt,'XT',[],[]);
	
	tstr = optCheck(figopt,'tstr',[],[]);
	skipStats = optCheck(figopt,'skipStats',false,[true false]);
	
    if plotShuf
        if skipStats
            nC = 2;
        else
            nC = 4;
        end
    else
        if skipStats
            nC = 1;
        else
            nC = 3;
        end
    end
    gArr = gobjects(nC,1);	
	
	f = specStruct.f;
    if ~isempty(fld1)
        data = specStruct.(fld1);
        fld1str = fld1;     
    else
       fld1str = 'orig';
       data = specStruct;
    end

    legstrs = cell(size(gArr));
    legstrs{1} = Lmap(fld1str);
    if ~skipStats
        legstrs{2} = Lmap('thr');
        legstrs{3} = Lmap('sig');
    end
    if plotShuf && ~skipStats
        legstrs{4} = Lmap('shuf');
    elseif plotShuf
        legstrs{2} = Lmap('shuf');
    end
    
    pow = data.(fld2);
    stats = data.stats.(fld2);
    k = 1;	
	
	if isempty(inputax)
        gArr(1) = plot(f,pow,'Color',clrsMap(fld1str));
    else
        gArr(1) = plot(inputax,f,pow,'Color',clrsMap(fld1str));
    end
    if ~isempty(lw)
        gArr(1).LineWidth = lw;
    end %if

	if dogrid; grid on; end
	
	hold on;

    if ~skipStats
        gArr(2) = plot(f,ones(size(f))*stats.sig_thresh_1T,'LineStyle','--','Color',clrsMap('thr'),'LineWidth',Wmap('thr'));
    
        if isempty(tstr)
            title([fld1str,' ',fld2])
        else
            title(tstr)
        end

        if ~isempty(lw)
            gArr(2).LineWidth = lw;
        end %if
    
        sigvec = stats.sigvec1TL;
	    hold on;
	    gArr(3) = plot(f,sigvec,'.','MarkerSize',mksz,'Color',clrsMap('sig'));
        if ~isempty(lw)
            gArr(3).LineWidth = lw;
        end %if
    end %if
	
    if plotShuf
        hold on;
        powS = data.comp.(['mnSh',fld2]);
        gArr(numel(gArr)) = plot(f,powS,'Color',clrsMap('shuf'));
        if ~isempty(lw)
            gArr(numel(gArr)).LineWidth = lw;
		end %if
    end	
	
	xlabel(flbl)
	ylabel('PSD')
	if isempty(XT)
        xticks([0:XTdt:XL(2)])
    else
        xticks(XT)
    end
	xlim(XL)
    if ~isempty(YL)
        ylim(YL)
    end	
    XTL = xticklabels();
    ax = gca;
    XTL(mod(cellnummat(XTL),XTLmod)~=0) = {''};
    if isempty(XT)
        ax.XTickLabels = XTL;  
    end	
    ax.TickDir = tickdir;
    ax.Box = dobox;
    xtickangle(0)
    if doleg
       legend(gArr,legstrs,'Location',legloc,'Box','off')
    end 

    if ~isempty(arrowfreq)
        yL = ax.YLim;  
        t = text(arrowfreq/XL(2),0,'\uparrow','FontSize',16,...
            'HorizontalAlignment','center','VerticalAlignment','bottom',...
            'Units','normalized','FontWeight','bold');  
    end
	

end %fn
