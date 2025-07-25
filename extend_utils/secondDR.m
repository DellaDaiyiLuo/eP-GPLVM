%% Dimension reduction on training data latents
% Include speed into consideration

tgrid = setopt.tgrid(2:400);
d_tc = [0 find(diff(tgrid')>1) numel(tgrid)];
xplot = result_la.xxsamp(2:400,:); % xinit; % xppca;%setopt.xplds;
speed = [];
for i=1:numel(d_tc)-1
    x_tmp = xplot(d_tc(i)+1:d_tc(i+1),:);
    x_sm = smoothdata(x_tmp,1,'gaussian',6);
    c = diff(x_sm);
    speed = [speed;c;c(end,:)];
end
figure;quiver3(xplot(:,1),xplot(:,2),xplot(:,3),speed(:,1),speed(:,2),speed(:,3))

speed_n = speed./std(speed,[],1);
xplot = [result_la.xxsamp(2:400,:) speed_n];

% ------
D = pdist(xplot);
Z = squareform(D);
Y = mdscale(Z,3);
%-------

k = 6;
display(['Number of neighbors to be searched: ' num2str(k)])
Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');

EdgeTable = zeros(size(Idx,1)*(k-1),2);
for i = 2:size(Idx,2)
    EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
end

H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
figure;plot(H)
L = laplacian(H);
[V,D] = eigs(L,5,'sa');

xx_dr2 = V(:,find(diag(D)>1e-10,3));
% figure;scatter(xx_dr2(:,1),xx_dr2(:,2),3,xx)
figure;scatter3(xx_dr2(:,1),xx_dr2(:,2),xx_dr2(:,3),3,xx)
