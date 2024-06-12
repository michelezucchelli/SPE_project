function ci = bootstrap_ci(data, statfun, gamma, nboot)
    % preallocate array for bootstrap estimates
    bootstat = zeros(nboot,1); 
    for r = 1:nboot
        % draw n elements with replacement from data
        bootSample = datasample(data,length(data)); 
        % compute desired statistic as a function of the replicates
        bootstat(r) = statfun(bootSample); 
    end
    % order the computed statistics
    bootstat = sort(bootstat); 
    % compute the confidence interval
    inf = floor(nboot * (1-gamma) / 2) + 1;
    sup = nboot + 1 - inf;
    ci = [bootstat(inf) bootstat(sup)];
end