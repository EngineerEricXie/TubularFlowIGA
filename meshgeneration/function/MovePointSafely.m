function [points,accepted_alpha,minimum_scaled] = MovePointSafely(points,elements,point_index,candidate,required_scaled)
% Backtrack a point move until every incident hex remains valid.
if point_index<1 || point_index>size(points,1)
    error('Point index %d is out of range',point_index);
end
if any(~isfinite(candidate))
    error('Candidate position for point %d is not finite',point_index);
end
incident=find(any(elements==point_index-1,2));
if isempty(incident)
    error('Point %d is not connected to any volume element',point_index);
end

original=points(point_index,:);
alphas=[1,0.5,0.25,0.125,0.0625,0.03125,0.015625,0];
accepted_alpha=NaN;
minimum_scaled=-inf;
for alpha=alphas
    points(point_index,:)=original+alpha*(candidate-original);
    [minimum_det,trial_scaled,bad_elements]=HexMeshQuality(points,elements,incident);
    if bad_elements==0 && minimum_det>0 && trial_scaled>=required_scaled
        accepted_alpha=alpha;
        minimum_scaled=trial_scaled;
        return;
    end
end

points(point_index,:)=original;
error('No valid position found for point %d; the unmodified incident mesh is invalid',point_index);
end
