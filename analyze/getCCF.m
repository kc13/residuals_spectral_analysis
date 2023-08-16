function [xc,lags,xcraw] = getCCF(src_delta,targ_delta,maxlag)
% getCCF.m - for normalized and raw CCF

assert(isequal(size(src_delta),size(targ_delta)) && iscell(src_delta) == iscell(targ_delta),...
	'mismatching dims and/or array type for src_delta and targ_delta')

if isvector(src_delta) && ~iscell(src_delta) % single continuous spike train
	% note xcorr convention of target in first arg
	[xcraw,lags] = xcorr(targ_delta,src_delta,maxlag,'none'); 
	xc = xcraw./((sum(src_delta) * sum(targ_delta))/numel(src_delta));
else  % trial per row
	nTr = size(src_delta,1);
	if ~iscell(src_delta) % num matrix
		xc_tr = cell2mat(arrayfun(@(t) xcorr(targ_delta(t,:),src_delta(t,:),maxlag,'none'),...
			[1:nTr]','UniformOutput',false));
		denom = (sum(src_delta(:)) * sum(targ_delta(:)))/numel(src_delta);
	else % cell matrix
		xc_tr = cell2mat(arrayfun(@(t) xcorr(targ_delta{t},src_delta{t},maxlag,'none'),...
			[1:nTr]','UniformOutput',false));
		denom = (sum([src_delta{:}]) * sum([targ_delta{:}]))/sum(cellnumel(src_delta)); %uses helper fn
	end %if
	xcraw = sum(xc_tr,1);
	xc = xcraw./denom;
	lags = [-maxlag:maxlag];	
end %if

end %fn