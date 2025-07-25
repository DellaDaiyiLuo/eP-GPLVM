load('track1_100ms.mat')
spikes = spikes(:,idx);
position = position(idx);
spikes = spikes(:,1000:end);
position = position(1000:end);
[pks,locs] = findpeaks(position,'MinPeakHeight',200,'MinPeakDistance',80);
locs_st = [locs,numel(position)];
[pks,locs] = findpeaks(-position,'MinPeakHeight',-50,'MinPeakDistance',80);
locs_end = locs;
fwd = [];
bwd = [];
for i=1:numel(locs_st)-1
    fwd = [fwd,locs_st(i):locs_end(i)];
    bwd = [bwd,locs_end(i):locs_st(i+1)];
end
figure;hold on;plot(position)
scatter(fwd,position(fwd))
scatter(bwd,position(bwd))

%-----------repeat for fwd and bwd -----------------%
pos_dir = position(fwd);
spk_dir = spikes(:,fwd);

pp = round(pos_dir/2);
tbl = tabulate(pp);
tc = zeros(size(spk_dir,1),size(tbl,1));
for i=1:size(tbl,1)
    if tbl(i,2)>0
        tc(:,i)=mean(spk_dir(:,pp==tbl(i,1)),2);
    end
end

tc = tc(:,tbl(:,2)>0);
tbl_ = tbl(tbl(:,2)>0,:);

% plot
[M, I] = max(tc,[],2);
[~, I_cell] = sort(I);
tc_n = tc./repmat(M,1,size(tc,2));
imshow(tc_n(I_cell,:))


tc_sm = smoothdata(tc,2,'gaussian',10);
figure;plot(tc(6,:));
hold on;plot(tc_sm(6,:))
title({'tc smoothed with Gaussian of window 10'})
R=corrcoef(tc_sm);
figure;image(R,'CDataMapping','scaled')

%------------------------%
tc_sm_fwd = tc_sm;
tbl_fwd = tbl_;

tc_sm_bwd = tc_sm;
tbl_bwd = tbl_;

%------------------------%
tc_sm = [tc_sm_bwd(:,end:-1:1),tc_sm_fwd];
tbl = [-1*tbl_bwd(end:-1:1,:);tbl_fwd];

load('track1_PBEs_4ms.mat')

nt = size(spikes,2);
spikes = double(spikes);
scaler = 1/25;
tc_sc = tc_sm.*scaler; %tc = tc_sm.*scaler_cell+0.0001;
tc_sc = tc_sc+min(nonzeros(tc_sc))/10;

% tc_scc = tc.*scaler_ratio+0.00001;
% tc_sc = smoothdata(tc_scc,2,'gaussian',10);

loglikelihood = -repmat(sum(tc_sc',2)',nt,1) + spikes'*log(tc_sc);
[Ls, xinitidx] = max(loglikelihood,[],2);
xinit_ratio = tbl(xinitidx,1);

figure
plot(xinit_ratio,'.')
xlim([4739,4793]);
xlim([4296,4350]);
%-------
save('pbe_bayes_llh.mat','loglikelihood','tbl','tc_sc');
matrix = exp(loglikelihood(4739:4793,:)');
% matrix = exp(loglikelihood(4739:4793,:)'); 
matrix_n = zeros(size(matrix));
for i=1:size(matrix,2)
    matrix_n(:,i) = matrix(:,i)/sum(matrix(:,i));
%     matrix_n(:,i) = (matrix(:,i)-min(matrix(:,i)))/(max(matrix(:,i))-min(matrix(:,i)));
end
figure;image(1:size(matrix,2),tbl(:,1)*2,matrix_n,'CDataMapping','scaled')
c = gray;
c = flipud(c);
colormap(c);
set(gca,'YDir','normal')
xlabel('time bin')
ylabel('Inbound      Outbound')
cb = colorbar;
cb.Label.String='Probability';
title('PBE example2')

%----------
tbl_idx = [size(tbl,1):-1:size(tbl_fwd,1)+1,1:size(tbl_fwd,1)];
figure;plot(tbl(tbl_idx,1))