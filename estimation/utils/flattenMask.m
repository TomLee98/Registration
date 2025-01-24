function mask = flattenMask(mask)
%FLATTENMASK This function flatten mask as 'F' order and transform as
%sparse matrix

arguments
    mask    (:,:,:)     {mustBeNonsparse}
end

d = unique(mask);
d(d==0) = [];       % remove background
d(isnan(d)) = [];   % remove possible uncentain components

row_mask = [];
col_mask = [];

for k = 1:numel(d)
    vd_loc = find(mask == d(k));   % F order, linear sub index
    row_mask = [row_mask; vd_loc]; %#ok<AGROW>
    col_mask = [col_mask; k*ones(size(vd_loc))]; %#ok<AGROW>
end
vol_mask = ones(size(row_mask));

% generate sparse matrix
mask = sparse(row_mask, col_mask, vol_mask, numel(mask), numel(d));
end

