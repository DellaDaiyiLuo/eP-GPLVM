function llh = comp_LLHtc(X, F, Xg, Fg, result, Y, tgrid)
%% Joint probability P(X,F,Y|t,Ftc, Xtc)

fftc = get_tc(Xg,Fg,X,result.rhoff,result.lenff);
rhoxx = result.rhoxx;
lenxx = result.lenxx;
rhoff = result.rhoff;
lenff = result.lenff;
sigma = result.sigma;
[nt,nneur] = size(F);

% P(X|t)
[Bfun, ~, ~, ~, ~, Kt] = prior_kernel(rhoxx,lenxx,nt,1,tgrid); 
uu = Bfun(X,1);
logKttrm = sum(log(diag(Kt)));
logp_x = -.5*logKttrm-.5*nt*log(2*pi)-.5*trace(uu'*uu);

% P(F|X)
fq = F-fftc;
logp_f_x = - .5*trace(fq'*fq)/sigma;

% P(Y|F)
logp_y_f = sum(-exp(F)+Y.*F-log(factorial(Y)),'all');

llh = logp_x+logp_f_x+logp_y_f;