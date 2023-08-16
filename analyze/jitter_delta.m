function [delta_jtr] = jitter_delta(spk_delta,jtr_int)
% jitter_delta.m - permute binary vector data within intervals 
% as in Amarasingham et al. (2012)
% 7/16/22 KC

%inputs: 
% spk_delta = binary vector, NxM matrix, or Nx1 cell array
% -> if matrix: N = # trials, M = # msec/trial
% -> 1 = spike occurred in the corresponding 1 msec bin
% -> if cell array: the Nth cell contains a vector of 0s/1s 
% 	 for Nth trial; allows for variable trial durations
% jtr_int = duration of jitter interval (msec)
% -> positive scalar
% -> requirement: length of spike delta vector or trial rows 
% must be evenly divisible by jtr_int

%outputs:
% delta_jtr = vector/matrix/cell array matching dims of spk_delta
% 0/1 entries randomly permuted within blocks of jtr_int msec

if isvector(spk_delta) && ~iscell(spk_delta) % single continuous spike train
	len = length(spk_delta);
	assert(mod(len,jtr_int)==0,'length % jitter interval != 0')
	nInt = len/jtr_int;
	delta_vec = toCol(spk_delta);
	rs = reshape(delta_vec,jtr_int,nInt);
	delta_jtr = reshape(cell2mat(arrayfun(@(x) randsample(rs(:,x),jtr_int),1:nInt, 'UniformOutput', false)),size(spk_delta));
else % trial per row
	if ~iscell(spk_delta) % num matrix
		[nTr,trlen] = size(spk_delta);
		assert(mod(trlen,jtr_int) == 0,'trial length % jitter interval != 0')
		nIntTr = trlen/jtr_int;
		rs = reshape(spk_delta',jtr_int,nIntTr*nTr);
		delta_jtr = reshape(cell2mat(arrayfun(@(x) randsample(rs(:,x),jtr_int),1:nIntTr*nTr,'UniformOutput',false)),trlen,nTr)';
	else % cell
		nTr = size(spk_delta,1);
		delta_jtr = cell(nTr,1);
		for tr = 1:nTr
			delta_vec = spk_delta{tr}';
			len = length(delta_vec);
			assert(mod(len,jtr_int)==0,'length % jitter interval != 0')
			nInt = len/jtr_int;
			rs = reshape(delta_vec,jtr_int,nInt);
			delta_jtr{tr} = reshape(cell2mat(arrayfun(@(x) randsample(rs(:,x),jtr_int),1:nInt, 'UniformOutput', false)),1,len);
		end %for
	end % if
end %if

end %fn