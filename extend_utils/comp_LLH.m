function llh = comp_LLH(result_la, yy, tgrid)
X = result_la.xxsamp;
F = result_la.ffmat;
Y = yy;
rhoxx = result_la.rhoxx;
lenxx = result_la.lenxx;
rhoff = result_la.rhoff;
lenff = result_la.lenff;
sigma = result_la.sigma;
[nt,nneur] = size(F);

% P(X|t)
[Bfun, ~, ~, ~, ~, Kt] = prior_kernel(rhoxx,lenxx,nt,1,tgrid); 
uu = Bfun(X,1);
logKttrm = sum(log(diag(Kt)));
logp_x = -.5*logKttrm-.5*nt*log(2*pi)-.5*trace(uu'*uu);

% P(F|X)
covfun = covariance_fun(rhoff,lenff,2);
Kff = covfun(X,X)+sigma*eye(size(X,1));
L = chol(Kff, 'lower');
logKftrm = sum(log(diag(L)));
logp_f_x = 0;
for n = 1:nneur
    squaretrm = F(:,n)'*(L'\(L\F(:,n)));
    logp_f_x = logp_f_x - .5*squaretrm - logKftrm - .5*nt*log(2*pi);
end

% P(Y|F)
logp_y_f = sum(-exp(F)+Y.*F-log(factorial(Y)),'all');

llh = logp_x+logp_f_x+logp_y_f;
