function r = GenROI(nid, center, radius)

arguments
    nid     (:,1)   cell
    center  (:,1)   cell
    radius  (:,1)   cell
end

if (numel(nid) ~= numel(center)) || (numel(center)~=numel(radius))
    throw(MException("CR2ROI:sizeNotMatch", ...
        "Input arguments size not match."));
end

N = max(cell2mat(cellfun(@(x)max(x), nid, "UniformOutput",false)));

r = cell(N, 1);

for zidx = 1:numel(nid)
    id_z = nid{zidx};
    for nidz = 1:numel(id_z)    % empty omitted
        if ~isnan(id_z(nidz)) && id_z(nidz)>=1   % certainty
            dr = [center{zidx}(nidz,:), radius{zidx}(nidz), zidx];
            r{id_z(nidz)} = [r{id_z(nidz)}; dr];
        end
    end
end

end

