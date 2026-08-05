function D = PointAlign(A,B,C)
% Project C onto the line through A and B.
direction=A-B;
length_squared=dot(direction,direction);
if ~isfinite(length_squared) || length_squared<=eps
    error('PointAlign requires two distinct finite line points');
end

k=dot(direction,C-B)/length_squared;
D=B+k*direction;
end
