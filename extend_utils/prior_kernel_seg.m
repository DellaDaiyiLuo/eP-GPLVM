function [BBwfun_seg, BBwTfun_seg, nu] = prior_kernel_seg(rhoxx,lenxx,nt,latentTYPE,tgrid, dt)
% prior_kernel_seg compute prior_kernel for each continuous segments
% Not speeding up the process. Jokes on me sorey!

% find continuous segments in tgrid
nu = nt;
d = [0,find(diff(tgrid')>dt*5),numel(tgrid')]; % if two continuous segments are 5 timebins apart we still construct the K_t together
BBwfun_array = cell(numel(d)-1,1);
BBwTfun_array = cell(numel(d)-1,1);

% compute prior_kernel for each segments
for i_seg = 1:numel(d)-1
    [BBwfun, BBwTfun, ~] = prior_kernel(rhoxx,lenxx,d(i_seg+1)-d(i_seg),latentTYPE,tgrid(d(i_seg)+1:d(i_seg+1)));
    BBwfun_array{i_seg} = BBwfun;
    BBwTfun_array{i_seg} = BBwTfun;
end

% compute function value by distributing into function arrays
BBwfun_seg = @(xx,invflag) dist_fun_array(xx,invflag,BBwfun_array,d);
BBwTfun_seg = @(xx,invflag) dist_fun_array(xx,invflag,BBwTfun_array,d);

function uu = dist_fun_array(xx,invflag,fun_array,d)
uu = zeros(size(xx));
for i_seg = 1:numel(d)-1
    part = d(i_seg)+1:d(i_seg+1);
    uu(part,:) = fun_array{i_seg}(xx(d(i_seg)+1:d(i_seg+1),:),invflag);
end
