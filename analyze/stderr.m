function [s] = stderr(data, varargin)
%stderr.m
% history
% 9/28/20 updated to add omitnan
% 2nd optional arg is dim
% 3rd option arg is nanflag
% 10/23/20: noticed header line said se, changed to stderr
% to match file name
% 5/7/23 updating so default dim is largest input dim

dropnan = false;
    if nargin >= 2 && ~isempty(varargin{1})
        dim = varargin{1};
        if ~isnumeric(dim)
            error('unrecognized second argument');
        end
        if nargin == 3
            if strcmp(varargin{2},'omitnan')
                dropnan = true;
            else
                error('unrecognized third argument');
            end
        end
    else
        [maxsz,maxI] = size(data);
        dim = maxI;
    end

   
    % note: using sample std -- this is different than what DK used
    if dropnan
        s = std(data,0,dim,'omitnan')./sqrt(sum(~isnan(data),dim));
    else
        s = std(data,0,dim)./sqrt(size(data,1));
    end
end

