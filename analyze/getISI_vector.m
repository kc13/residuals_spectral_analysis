function [ISI] = getISI_vector(delta_vec)

	ISI = diff(find(delta_vec));

end %fn