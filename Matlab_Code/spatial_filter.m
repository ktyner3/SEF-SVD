function [output] = spatial_filter(data, points, distance)
% data = data you wish to spatially filter
% points = 3D source locations
% distance = spatial parameter for filter
Mdl = KDTreeSearcher(points);
Idx = rangesearch(Mdl,points,distance);
output = cell(length(Idx),1);
for j = 1:length(Idx)
    pts = Idx{j,1};
    output{j,1} = data(pts);
    output{j,1} = mean(output{j,1});
end

output = cell2mat(output);

end
