function B = RotateSurface(A,norm_old,norm_new)
% Apply the shortest rotation that maps norm_old onto norm_new.
old_length=norm(norm_old);
new_length=norm(norm_new);
if ~isfinite(old_length) || old_length<=eps || ~isfinite(new_length) || new_length<=eps
    error('RotateSurface requires finite, nonzero normals');
end
norm_old=norm_old/old_length;
norm_new=norm_new/new_length;
cosine=max(-1,min(1,dot(norm_old,norm_new)));
if cosine>=1-1.0e-14
    B=A;
    return;
end
if cosine<=-1+1.0e-14
    [~,least_aligned]=min(abs(norm_old));
    helper=zeros(1,3);
    helper(least_aligned)=1;
    rotation_axis=cross(norm_old,helper);
    alpha=pi;
else
    rotation_axis=cross(norm_old,norm_new);
    alpha=acos(cosine);
end
B=RotateAroundAxis(A,rotation_axis,alpha);
end
