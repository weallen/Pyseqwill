mm9_chr_lengths =  [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296, 15902555];
mm9_chr_bnd = floor(mm9_chr_lengths / 200);

DATA_H5 = '~/Documents/Neuroscience/barnea_lab/rna_seq_experiment/stavros_chip/all_data.h5';
omp_hmedip = hdf5read(DATA_H5, '/seqdata/omp_hmedip');
chr1_omp_hmedip = omp_hmedip(1:mm9_chr_bnd(1));

chr1_hme_model = hmmFit(chr1_omp_hmedip, 2, 'discrete');
