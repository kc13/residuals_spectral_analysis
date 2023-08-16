function [sim_spk_t,p,p_raw] = simulate_spk(len, prob_spk_base, nr,k, f_osc, p_osc, osc_offset)
%simulate_spk.m 
% --> simulates spike times for a provided interval
% --> following from Rivlin-Etzion et al., 2006, JNP

sim_spk_t = nan(len,1);

% check for osc_offset
if ~(exist('osc_offset','var'))
	osc_offset = 0;
end

% expand scalar args to vector
p_osc = expand_param(p_osc,len);
prob_spk_base = expand_param(prob_spk_base,len);
nr = expand_param(nr,len);
k = expand_param(k,len);

% init n to be in last bin of RP, loop will start by incrementing by 1
n = nr(1);

p = nan(len,1);
p_raw = nan(len,1);
% loop over ms time bins
for t = 1:len

	n = n + 1; % updating timesteps since last spike, assuming present bin should be counted
	    
	% probability considering base and oscil. 
	p_raw(t) = prob_spk_base(t) + p_osc(t)*sin(2*pi*f_osc*(t+osc_offset)*0.001);
	
    if n > nr(t)
        p(t) = p_raw(t); % outside of R.P.
    else     
        p(t) = k(t)^(nr(t) + 1 - n)*p_raw(t); % considers RP
    end   
	
	assert(p(t) >= 0 && p(t) <= 1,'invalid simulation probabilty');
	
    % runs 1 Bernoulli trial
    sim_spk_t(t) = binornd(1,p(t)); 

	% reset n if spike
	if sim_spk_t(t) == 1
		n = 0;
	end
	
end % t 


    function [xe] = expand_param(x,len)
        numE = numel(x);
        if numE == 1
            xe = repmat(x,len,1);
        else
            assert(numE == len,'number of elements in a vector does not match requested length')
            xe = x;
        end %if
    end %fn

end  % fn 