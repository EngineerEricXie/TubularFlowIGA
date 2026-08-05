function Q = Projection(A,B,C,P)
% Project P onto the plane defined by A, B, and C.
normal=cross(B-A,C-A);
normal_length=norm(normal);
if ~isfinite(normal_length) || normal_length<=eps
    error('Projection requires three finite, non-collinear plane points');
end
normal=normal/normal_length;

Q=P-dot(P-A,normal)*normal;
end
