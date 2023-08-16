function [resmat,mdls,rpInfo] = getLastSpkResiduals(data,resopt)
%getLastSpkResiduals.m

concatWins = optCheck(resopt,'concatWins',false,[true false]);
oneModel = optCheck(resopt,'oneModel',true,[true false]);
model_order = optCheck(resopt,'order','auto',[]);

if isnumeric(model_order)
    rpInfo.RPend = model_order;
end

runSplines = optCheck(resopt,'runSplines',false,[true false]);
res_options = optCheck(resopt,'res_options',statset('fitglm'),[]);

distname = optCheck(resopt,'distname','poisson',[]);
DF = optCheck(resopt,'DispersionFlag',false,[true false]);

if runSplines 
   splinePrefs = resopt.splinePrefs; 
end

run_parfor = optCheck(resopt,'run_parfor','true',[true false]);
verbose = optCheck(resopt,'verbose','true',[true false]);

if ~verbose
    % temp disable
    warning('off','all')
end

assert(~(concatWins && ~oneModel),'this concatWins / oneModel combo not yet supported')

[nW,len] = size(data);

if strcmp(model_order,'auto') % run RP estimation
	if ~concatWins % discontiguous data
		if ~oneModel % separate RPs estimated per window
			order = nan(nW,1);
			p = cell(order,1);
			for w = 1:nW
			   rpInfoW = runRPfinder_vector(data(w,:));
			   order(w) = rpInfoW.orderS;
			   p{w} = rpWinfoW.p;
			end			
			order_vec = order; 
			rpInfo.orderS = order;
			rpInfo.p = p;
		else % one RP
			[rpInfo] = runRPfinder(data);
			order = rpInfo.RPend;
			order_vec = repmat(order,nW,1);
		end
	else % single unbroken vector
		data_vec = reshape(data',[],1);
		[rpInfo] = runRPfinder_vector(data_vec);  
		order = rpInfo.RPend;
	end %if ~concat wins
else % hard coded order
    order = model_order;
    rpInfo.order = order;	
end %if strcmp auto

if runSplines && order >= splinePrefs.minSplPts
    splinePrefs.lag = order;
    splinePrefs.lastknot = splinePrefs.lag; 
	splinePrefs.c_pt = getCPT(splinePrefs.lag,splinePrefs.method,splinePrefs.brkpt);
	rpInfo.splinePrefs = splinePrefs; 
    Spl = getSplines(splinePrefs.method,splinePrefs.lag,splinePrefs.c_pt,splinePrefs.s);
	confirmSplines = true;
else
	confirmSplines = false;
end

% preallocate
if ~concatWins
    resmat = nan(size(data));
    if oneModel
        xpredAll = cell(nW,1);
        yvecAll = cell(nW,1);
    end
end

% this step builds and runs the models
if ~concatWins
    for w = 1:nW
        ix = order_vec(w)+1:len; 
        ymat = data(:,ix);
        yvec = ymat(w,:)';
        dW = data(w,:);
        dWcAll = arrayfun(@(x) circshift(dW,x,2),[1:order_vec(w)],'UniformOutput',false);
        xmatAllRaw = cell2mat(cellfun(@(x) x(:,ix)', dWcAll,'UniformOutput',false));
        xmatAllTag = arrayfun(@(x) tagFirst(xmatAllRaw(x,:)),[1:size(xmatAllRaw,1)]','UniformOutput',false); 
        xpred = cell2mat(xmatAllTag);
        if ~oneModel
            % splines only implemented for the single RP model
			assert(~confirmSplines,'run splines not yet ready for this case')
            mdls{w} = fitglm(xpred,yvec,'linear','Distribution',distname,'Options',res_options,'DispersionFlag',DF);
            resmat(w,ix) = mdls{w}.Residuals.Raw; 
        else
            xpredAll{w} = xpred;
            yvecAll{w} = yvec;
        end %if ~oneModel
    end %w    
else % can work with vector
    datacol = reshape(data',[],1);
    ycol = datacol(order+1:end); 
    xcols_raw = cell2mat(arrayfun(@(x)...
        toCol(datacol(order-x+1:end-x)),[1:order],'UniformOutput',false));
    xcols = double(cell2mat(arrayfun(@(x) tagFirst(xcols_raw(x,:)), ...
        [1:size(xcols_raw,1)]','UniformOutput',false)));
end %if ~concatWins

if oneModel
    if ~concatWins
		assert(~confirmSplines,'run splines not yet ready for this case')
        mdls = fitglm(vertcat(xpredAll{:}),vertcat(yvecAll{:}),'linear','Distribution',distname,'Options',res_options,'DispersionFlag',DF);
        resmat(:,order+1:end) = reshape(mdls.Residuals.Raw,size(data,2)-order,nW)';          
    else
		if confirmSplines
            scols = xcols*Spl;
            disp('running fitglm')
			mdls = fitglm(scols,ycol,'linear','Distribution',distname,'Options',res_options,'DispersionFlag',DF); 
		else
			disp('running fitglm')
			mdls = fitglm(xcols,ycol,'linear','Distribution',distname,'Options',res_options,'DispersionFlag',DF);
		end
        resmat = reshape(vertcat(nan(order,1),mdls.Residuals.Raw),size(data,2),nW)';
    end %if
end % if

% restore warnings
if ~verbose
    warning('on','all')
end


end %fn
