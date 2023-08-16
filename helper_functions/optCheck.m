function [val] = optCheck(opts,fld,default,validSet)
%optCheck.m

if ~isfield(opts,fld)
    val = default;
elseif isempty(validSet) || ismember(opts.(fld),validSet)
    val = opts.(fld);
else    
    error('unrecognized value %s for field %s',string(opts.(fld)),fld)
end %if

end %fn

