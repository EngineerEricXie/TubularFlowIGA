function P = RotateAroundAxis(Q,axis,angle)
% Rotate points Q around axis by angle radians.
axis_length=norm(axis);
if abs(angle)<=1.0e-14
    P=Q;
    return;
end
if ~isfinite(axis_length) || axis_length<=eps
    error('RotateAroundAxis requires a finite, nonzero axis');
end
axis=axis/axis_length;
K=[0 -axis(3) axis(2);
   axis(3) 0 -axis(1);
   -axis(2) axis(1) 0];
I=eye(3);
R=I+sin(angle)*K+(1-cos(angle))*K*K;
P=Q*R';
end
