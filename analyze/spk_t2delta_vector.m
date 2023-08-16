function [spk_delta] = spk_t2delta_vector(spk_t,data_len)
%spk_t2delta_vector.m
%5/22/22

spk_delta = zeros(data_len,1); % column to save memory

% this is spk_t2delta convention
tmp_spk = 1+round(1000*spk_t);

% will cut of anything trailing if exceeds bounds after rounding
spk_delta(tmp_spk(tmp_spk <= data_len)) = 1;

rmv = numel(spk_t) - nnz(spk_delta);

if rmv > 0 
	warning('%d trailing spikes removed in conversion to delta\n',rmv)
end

end %fn