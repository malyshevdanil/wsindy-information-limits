function Cfs = phi_weights(phifun,m,maxd)
  % Exact derivatives of the author default exp(-9/(1-x^2)) bump.
  assert(maxd<=1);
  x=linspace(-1,1,2*m+1); x=x(2:end-1);
  f=phifun(x);
  assert(max(abs(f-exp(-9./(1-x.^2))))<1e-14);
  Cfs=zeros(maxd+1,2*m+1);
  Cfs(1,2:end-1)=f;
  if maxd==1
    Cfs(2,2:end-1)=f.*(-18*x./(1-x.^2).^2);
  end
  Cfs=Cfs/norm(Cfs(1,:),2);
end
