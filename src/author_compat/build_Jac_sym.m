function Jac_mat = build_Jac_sym(features,xobs)
  % Exact polynomial differentiation replaces MATLAB symbolic plumbing only.
  global benchmark_feature_powers;
  [M,d] = size(xobs);
  J = sum(cellfun(@length,features));
  assert(size(benchmark_feature_powers,1)==J);
  Jac_mat = zeros(J,d,M);
  for j=1:J
    p = benchmark_feature_powers(j,:);
    for r=1:d
      if p(r)>0
        q=p; q(r)=q(r)-1;
        v=p(r)*prod(xobs.^q,2);
        Jac_mat(j,r,:)=reshape(v,1,1,M);
      end
    end
  end
end
