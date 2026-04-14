function SWI = swi(Wthresh)
% Computes small-world index of Wthresh

L = weight_conversion(Wthresh, 'lengths');            % taking the inverse of matrix elements
D = distance_wei(L);                            % applying Dijkstra's algorithm
% compute characteristic path length and global efficiency
[lambdaW, ~] = charpath(D, 0, 0);
% random matrix
% WthreshR = threshold_proportional(Wr, threshold);  % Applying proportional threshold
WthreshR = rand(size(Wthresh));
Lr = weight_conversion(WthreshR, 'lengths');            % taking the inverse of matrix elements
Dr = distance_wei(Lr);                            % applying Dijkstra's algorithm
% compute characteristic path length and global efficiency
[lambdaWr, ~] = charpath(Dr, 0, 0);
% -------------- Clustering coefficient ---------------------
% Compute clustering coefficient
clust_coeff = nanmean(clustering_coef_bu(Wthresh));
clust_coeff_r = nanmean(clustering_coef_bu(WthreshR));
SWI = (clust_coeff / clust_coeff_r) / (lambdaW / lambdaWr);
end

% EOF