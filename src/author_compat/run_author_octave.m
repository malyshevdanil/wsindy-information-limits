function run_author_octave(input_path,output_path,repo_path,compat_path,method)
  warning('off','Octave:shadowed-function');
  pkg load statistics;
  pkg load signal;
  addpath(repo_path);
  data=load(input_path); xobs=data.y; tobs=data.t(:);
  started=tic;
  if strcmp(method,'wsindy')
    % Latest pinned complete author pipeline: FFT width, gradient placement,
    % automatic smoothing, scaling and MSTLS threshold selection.
    [coef,~,~,loss,~,elapsed,~,grids,pts,mts,Gs,bs,M_diag,RTs,Theta,tags,~,~,vs,filters,xsmooth] = ...
      wsindy_ode_fun(xobs,tobs,[],1:3,[],zeros(0,2),{},1,1,1e-16,-1,...
      1,1,1,2,0,10.^linspace(-4,0,100),0,.8,.8,[1;0],[1;0],0,2,ceil(length(tobs)/100));
    result=struct('coefficients',coef,'tags',tags,'mts',mts,'pts',pts,...
      'grid_counts',cellfun(@length,grids),'filter_lengths',cellfun(@length,filters),...
      'grid_indices',{grids},'loss',{loss},'elapsed',elapsed);
  elseif strcmp(method,'wendy')
    addpath(compat_path);
    global benchmark_feature_powers;
    benchmark_feature_powers=[0,1;1,0;0,1;2,1];
    features={{@(x1,x2)x2},{@(x1,x2)x1,@(x1,x2)x2,@(x1,x2)x1.^2.*x2}};
    wendy_snf_params;
    [M,nstates]=size(xobs); phifun=phifuns{1}; meth='mtmin';
    mt_params=2.^(0:3); K_max=5000; K_min=4;
    mt_max=max(floor((M-1)/2)-K_min,1);
    mt_min=rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
    mt_cell=cellfun(@(x,y)[x,{y}],repmat({{phifun,meth}},4,1),num2cell(mt_params(:)),'uni',0);
    % Dummy diagnostic vector: actual coefficients are withheld from fitting.
    true_vec=ones(4,1); use_true=0;
    [w,res,~,~,~,its,~,V,Vp,Theta,mt,~,Jac,G,b,RT,stdW,mseW,CovW] = ...
      wendy_fcn(xobs,tobs,features,true_vec,toggle_smooth,mt_cell,mt_min,...
      mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,w0,optim_params,...
      err_norm,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it,use_true);
    % Feature order is x2 in equation 1, then x1,x2,x1^2*x2 in equation 2.
    coef=zeros(9,2); tags=data.powers;
    for j=1:4
      r=find(all(tags==benchmark_feature_powers(j,:),2));
      eq=1+(j>1); coef(r,eq)=w(j);
    end
    result=struct('coefficients',coef,'tags',tags,'mts',mt,'mt_min',mt_min,...
      'basis_counts',cellfun(@(v)size(v,1),V),'iterations',size(its,2),...
      'iterates',its,'sigma_estimates',arrayfunvec(xobs,@estimate_sigma,1),...
      'covariance',CovW,'residual_norm',norm(res(:,end)));
    [~,pv,~]=swtest(res(:,end)); result.final_normality_pvalue=pv;
  else
    error('Unknown method');
  end
  result.total_seconds=toc(started);
  fid=fopen(output_path,'w'); fprintf(fid,'%s',jsonencode(result)); fclose(fid);
end
