% try to use empirical bayesian methods to determine threshold, with a 
% certain FDR

function enrich = genomeEnrichment(genome_coverage)
    empirical_mean = mean(genome_coverage);
    threshold = 0;
    for i=1:40
        if poisspdf(i, empirical_mean) < 1E-4
            threshold = i;
            break;
        end
    end  
    enrich = ones(1, length(genome_coverage));
    for i=1:length(genome_coverage)
        if genome_coverage(i) >= threshold
            enrich(i) = 2;
        end
    end
end
