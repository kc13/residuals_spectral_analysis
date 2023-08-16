function [T] = tagFirst(V)

T = zeros(size(V));
ix = find(V);
if ~isempty(ix)
    T(ix(1)) = V(ix(1));
end

end