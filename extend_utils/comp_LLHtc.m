function llh = comp_LLHtc(X, F, Xg, Fg, result, Y, tgrid,varargin)
%% Joint probability P(X,F,Y|t,Ftc, Xtc)

p = inputParser;
addParameter(p,'cell',0)
parse(p,varargin{:});

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

fq = F-fftc;

if p.Results.cell
    % P(F|X)
    logp_f_x = - .5*diag(fq'*fq)/sigma;
    % P(Y|F)
    logp_y_f = sum(-exp(F)+Y.*F-log(factorial(Y)),1);

    llh = [logp_f_x(:),logp_y_f(:)];
else
    % P(F|X)
    logp_f_x = - .5*trace(fq'*fq)/sigma;
    % P(Y|F)
    logp_y_f = sum(-exp(F)+Y.*F-log(factorial(Y)),'all');
    
    llh = logp_x+logp_f_x+logp_y_f;
end
