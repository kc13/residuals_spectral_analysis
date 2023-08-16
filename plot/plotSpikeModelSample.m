function [ax] = plotSpikeModelSample(data,msecs,category,varargin)
%plotPosterSpikeSample.m 2/21/23
% 6/19/23 updating to allow more formatting flexibility
% 7/2/23 updating to says "oscillation" in full    
    if ~isempty(varargin)
        figopt = varargin{1};
        add_annotate = optCheck(figopt,'add_annotate',false,[true false]);
        axftsz = optCheck(figopt,'axftsz',[],[]);
        lwInc = optCheck(figopt,'lwInc',0,[]);
        tftsz = optCheck(figopt,'tftsz',[],[])
    else
        add_annotate = false;
        axftsz = [];
        lwInc = 0;
        tftsz = [];
    end %if

    tickdir = 'out';
    dobox = 'off';

    xUB = max(msecs);
    xUBR = round(xUB,-2);
    xtint = 50;
    xt = [0:xtint:xUBR];

    clrsMap = containers.Map();
    clrsMap('spikes') = [0 0 1]; %[1 1 1]*0.6; %same as orig in plotSpikeSample
    clrsMap('oscil') = [0 0.65 0.45]; 
    clrsMap('rp') = [122 0 25]/255;
    clrsMap('pspk') = clrsMap('oscil').*clrsMap('rp');

    titleMap = containers.Map();
    titleMap('spikes') = 'spike train';
    titleMap('oscil') = 'steady state (base FR + oscillation)';
    titleMap('rp') = 'recovery period';
    titleMap('pspk') = 'spike probability = SS x RP';

    xlblMap = containers.Map();
    xlblMap('spikes') = 'data segment timestamp (1 ms bins)'; %'timestamp (1 ms bins)'; 
    xlblMap('oscil') = ''; % blank for now b/c redundant w/ spikes at end
    xlblMap('rp') = '';  % blank for now b/c redundant w/ spikes at end
    xlblMap('pspk') = '';

    ylblMap = containers.Map();
    ylblMap('spikes') = 'occurrence';
    ylblMap('oscil') = 'p_{SS}(t) (zoomed)';  % this will only make sense if show eq
    ylblMap('rp') = 'p_{RP}(t)'; 
    ylblMap('pspk') = 'p_{spk}(t) (zoomed)';

    ylimMap = containers.Map();
    ylimMap('spikes') = [0 1];
    ylimMap('oscil') = [0 0.3]; % from cosyne: [0 0.06]; 
    ylimMap('rp') = [0 1]; % will need to clarify
    ylimMap('pspk') = ylimMap('oscil');
    

    if isequal(category,'spikes')
        % all matching fig 2 here
        s = stem(msecs,data(msecs),'Color',clrsMap(category));
        s.Marker = 'o';
        s.ShowBaseLine = 'off';
        s.MarkerFaceColor = s.Color;
        s.MarkerEdgeColor = 'none';
        s.MarkerSize = 1;
        s.LineWidth = s.LineWidth+lwInc;
        ax = gca;
    else
        p = plot(msecs,data(msecs),'Color',clrsMap(category));
        p.LineWidth = p.LineWidth+lwInc;
        ax = gca;
    end %if

    if add_annotate
        switch category
            case 'oscil'
                disp('here')
                keyboard
        end %sw
    end %if    

    ax.YLim = ylimMap(category);
    ax = gca;
    ax.TickDir = tickdir;
    ax.Box = dobox;
    ax.XLim = [0 max(msecs)]; 
    xlabel(xlblMap(category))
    ylabel(ylblMap(category))
    t = title(titleMap(category));
    if ~isempty(tftsz)
        t.FontSize = tftsz;
    end
    if ~isempty(axftsz)
        ax.FontSize = axftsz;
    end
    %ax.FontSize = 11;
    
    ax.LineWidth = ax.LineWidth+1;
end
%%