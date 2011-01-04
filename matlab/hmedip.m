mm9_chr_lengths =  [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296, 15902555];
mm9_chr_bnd = ceil(mm9_chr_lengths / 200);

DATA_H5 = '~/experiment/experiment/stavros_data/all_data.h5';
omp_hmedip = hdf5read(DATA_H5, '/seqdata/omp_hmedip');
data = {genomeEnrichment(omp_hmedip)}; 
fprintf('Found enrichment\n');
omp_hme_model = hmmFit(data, 3, 'discrete');
fprintf('Fit model\n');
omp_hme_path = hmmMap(omp_hme_model, data{1});
fprintf('Found path\n');
fprintf('Found peaks\n');
windows = loadWindows(DATA_H5);
writePathToWigFile('omp_hme.wig', omp_hme_path, windows);
save('omp_hme_model','omp_hme_model');
save('omp_hme_path', 'omp_hme_path');
