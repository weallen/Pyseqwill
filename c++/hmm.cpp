#include "hmm.h"

void
load_seq_data() 
{
    std::string file_name("/gpfs/home/wallen/experiment/experiment/stavros_data/all_data.h5");
    H5::H5File file(file_name, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("/seqdata/windows");
    H5::DataSpace filespace = dataset.getSpace();
    hsize_t dims[2];
    int rank = filespace.getSimpleExtentNdims();
    int ndims = filespace.getSimpleExtentDims(dims, NULL);
    int nrows = dims[0];
    std::cout << "rank " << rank << ", dimensions " << 
        (unsigned long)(dims[0]) << " x " << (unsigned long)(dims[1]) << std::endl; 
}

void
HMM::init()
{
    boost::mt19937 rng(static_cast<unsigned int>(time(0)));
    boost::normal_distribution<float> noise(0, 1);
    boost::variate_generator<boost::mt19937, 
                             boost::normal_distribution<float> > gauss(rng, noise);
    for (unsigned i = 0; i < trans_probs_.size1(); ++i) {
        for (unsigned j = 0; j < trans_probs_.size2(); ++j) {
            trans_probs_(i, j) = gauss();
        }
    }
    for (unsigned i = 0; i < emission_probs_.size1(); ++i) {
        for (unsigned j = 0; j < emission_probs_.size2(); ++i) {
        }
    }
}
