#include <H5Cpp.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <string>

using namespace boost::numeric::ublas;

class WindowData {
public:
    WindowData(unsigned long nrows, int* data)
        : num_rows_(nrows), data_(data) {}

    unsigned long num_rows() const { return num_rows_; }
    int* data() { return data_; }

private:
    unsigned long num_rows_;   
    int* data_; 
};

class HMM {
public:
    HMM(int M, int K) 
        : M_(M), K_(K), trans_probs_(M_, M_), emission_probs_(M_, K_) 
        {
            init();
        }

    ~HMM() {}
    
    const matrix<float>& trans_probs() const { return trans_probs_; }
    const matrix<float>& emission_probs() const { return emission_probs_; }
    
private:
    void init();
    int M_;
    int K_;
    matrix<float> trans_probs_;
    matrix<float> emission_probs_; 
};

void load_seq_data();

