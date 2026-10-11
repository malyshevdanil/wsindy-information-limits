function C = pagemtimes(A,B)
  % Octave bridge for the two-input MATLAB call used by the pinned WENDy.
  % The benchmark has one batched operand; no transpose modes are used.
  n = max(size(A,3),size(B,3));
  C = zeros(size(A,1),size(B,2),n);
  for k=1:n
    C(:,:,k) = A(:,:,min(k,size(A,3))) * B(:,:,min(k,size(B,3)));
  end
end
