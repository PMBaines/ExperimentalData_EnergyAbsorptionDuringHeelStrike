function [mean_total,mean_min,mean_max,within_std,between_std,means] = function_GetStats(data,matrix,subjects)

D = (data.*matrix)';
D = D(:,subjects);

for index = 1:length(subjects)
    
    V = D(:,index);
    V = V(V~=0);
    V = V(isfinite(V));
    means(index) = mean(V);
    stds(index) = std(V);
end

mean_total = mean(means);
within_std = mean(stds);
between_std = std(means);
mean_min = min(means);
mean_max = max(means);

end