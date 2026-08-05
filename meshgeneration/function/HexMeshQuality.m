function [minimum_det,minimum_scaled,bad_elements,first_bad_element] = HexMeshQuality(points,elements,element_indices)
% Evaluate VTK hexahedra at their corners and four Gauss points per axis.
if nargin<3
    element_indices=1:size(elements,1);
end
samples=[-1,-0.861136311594053,-0.339981043584856,...
    0.339981043584856,0.861136311594053,1];
signs=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;...
    -1 -1 1;1 -1 1;1 1 1;-1 1 1];
minimum_det=inf;
minimum_scaled=inf;
bad_elements=0;
first_bad_element=0;

for element_index=element_indices(:)'
    node_indices=elements(element_index,:)+1;
    if any(node_indices<1) || any(node_indices>size(points,1))
        error('Element %d contains an out-of-range point index',element_index-1);
    end
    coordinates=points(node_indices,:);
    element_bad=false;
    for r=samples
        for s=samples
            for t=samples
                derivatives=zeros(8,3);
                derivatives(:,1)=0.125*signs(:,1).*(1+signs(:,2)*s).*(1+signs(:,3)*t);
                derivatives(:,2)=0.125*signs(:,2).*(1+signs(:,1)*r).*(1+signs(:,3)*t);
                derivatives(:,3)=0.125*signs(:,3).*(1+signs(:,1)*r).*(1+signs(:,2)*s);
                jacobian=coordinates'*derivatives;
                determinant=det(jacobian);
                scale=prod(vecnorm(jacobian,2,1));
                if ~isfinite(scale) || scale<=eps
                    scaled=-inf;
                else
                    scaled=determinant/scale;
                end
                minimum_det=min(minimum_det,determinant);
                minimum_scaled=min(minimum_scaled,scaled);
                element_bad=element_bad || ~isfinite(determinant) || determinant<=0;
            end
        end
    end
    if element_bad
        bad_elements=bad_elements+1;
        if first_bad_element==0
            first_bad_element=element_index;
        end
    end
end
end
