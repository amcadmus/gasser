#ifndef __wanghan_ToolBox_h__
#define __wanghan_ToolBox_h__
#include<vector>

typedef double value_type;

namespace RandomGenerator_MT19937{
    void init_by_array(unsigned long init_key[], int key_length);
    void init_genrand(unsigned long s);
    unsigned long genrand_int32(void);
    long genrand_int31(void);
    double genrand_real1(void); // in [0,1]
    double genrand_real2(void); // in [0,1)
    double genrand_real3(void); // in (0,1)
    double genrand_res53(void);
    void genrand_Gaussian (double mu, double sigma, double * rand1);
    void genrand_Gaussian (double mu, double sigma, double * rand1, double * rand2);
}




#endif
