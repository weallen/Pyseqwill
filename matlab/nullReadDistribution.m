% typical num_rand is 10
function expected = nullReadDistribution(num_reads, genome_len, window_size)
    rand_genome = zeros(1, genome_len);
    for i=1:num_reads
        win = ceil(randi(genome_len) / window_size);
        rand_genome(win) = rand_genome(win) + 1;
    end
    expected = mean(rand_genome);
end
