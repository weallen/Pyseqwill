function enrich = chromenrichment (chrom)
    empirical_mean = mean(chrom);
    for i=1:15
       if poisscdf(i, empirical_mean) < 1E-4
           threshold = i;
       end
    end
    enrich = zeros(length(chrom));
    for i=1:length(chrom)
        if chrom(i) > treshold;
            enrich(i) = 1;
        else 
            enrich(i) = 0;
        end
    end
end

