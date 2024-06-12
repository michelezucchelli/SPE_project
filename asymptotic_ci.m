function ci = asymptotic_ci(data, gamma)
    % the asymptotic formula cannot be applied if the number of sample is too small
    if length(data) < 30
        disp('The number of data in the dataset is too small')
        ci = [];
        return 
    end
    % (1+gamma)/2 -th quantile of the Normal distribution
    etha = norminv((1+gamma)/2);
    ci = mean(data) + [-etha etha]*std(data)/sqrt(length(data));
end