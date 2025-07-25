function [tbl_id, comp_id,k_log,n_debris_log] = find_components(xplot, ratio, plot_graph)
% Determine number of components in given data
%
% xplot: (#samples, #features), to be devided into disconnected components
% ratio: scalar, ratio between #(the smallest k to be search at the beginning) and #samples
% plot_graph: 1, plot graph

debris_size = 0.03;
k = max([3, floor(ratio*size(xplot,1))]);
nt = size(xplot,1);
n_debris = 1;
k_log = [];
n_debris_log = [];
n_comp_log = 10;
while n_comp_log(end)>1%n_debris>0 %
    k_log = [k_log k];
    display(['Number of neighbors to be searched: ' num2str(k)])
    Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');

    EdgeTable = zeros(size(Idx,1)*(k-1),2);
    for i = 2:k
        EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
    end

    H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
    comp_id = conncomp(H);
    tbl_id = tabulate(comp_id);
    n_debris = numel(find(tbl_id(:,2)<nt*debris_size)); % Any connected component should have more than 20 points
    n_debris_log = [n_debris_log n_debris];
    display(['number of debirs components: ' num2str(n_debris)])
    n_comp_log = [n_comp_log,size(tbl_id,1)];
    k = k+1;
end

n_comp_log = n_comp_log(2:end);
display(['Number of disconnected components: ' num2str(n_comp_log(end))])
n_comp = mode(n_comp_log);
k = k_log(find(n_comp_log==n_comp,1,'last'));
% k = k-2;
Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');
EdgeTable = zeros(size(Idx,1)*(k-1),2);
for i = 2:k
    EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
end
H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
comp_id = conncomp(H);
tbl_id = tabulate(comp_id);

if plot_graph==1
    figure;plot(H);
    title({'k-nearest neighbor graph' ['k=' num2str(k) ', n\_comp=' num2str(size(tbl_id,1))]})
    figure
    hold on
    plot(k_log,n_debris_log)
    plot(k_log,n_comp_log)
    xlabel('k')
    ylabel('number of components')
    legend({'residual','all'})
    title(['number of latent dimension: ' num2str(size(xplot,2))])
end

display(nt)
% while size(tbl_id,1)>1
%     k_log = [k_log k];
%     display(['Number of neighbors to be searched: ' num2str(k)])
%     Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');
% 
%     EdgeTable = zeros(size(Idx,1)*(k-1),2);
%     for i = 2:k
%         EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
%     end
% 
%     H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
%     comp_id = conncomp(H);
%     tbl_id = tabulate(comp_id);
%     n_debris = numel(find(tbl_id(:,2)<20)); % Any connected component should have more than 20 points
%     n_debris_log = [n_debris_log n_debris];
%     display(['number of debirs components: ' num2str(n_debris)])
%     k = k+1;
% end
end