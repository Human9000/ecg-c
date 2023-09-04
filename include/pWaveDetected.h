#ifndef PWAVEDETECTED_H
#define PWAVEDETECTED_H


#define type_0 0
#define type_1 1
#define type_2 2
#define type_3 3
#define type_4 4
#define type_5 5
#define type_6 6
//enum  pType{type_0,   type_1,   type_2,type_3,type_4,type_5,type_6};

// type_0    no
// type_1    positive
// type_2    inverse （negative）junctional
// type_3    double wave (positive negative) type 3
// type_4   double wave (negative positive)
// type5    ectopic P wave
// type_6   I AVB P wave

typedef struct {
    int onset;
    int offset;
    int type;
    float coeff;
    int Peak;
    float weight;
} pwav_data_t;


typedef struct {
    int onset;
    int offset;
    int type;
    int peak2;
    int peak1;
    float weight;
} pwav_data_t2;


//class pWaveDetected
//{


//public:

//    pWaveDetected();
//    ~pWaveDetected();


//public:

int pWaveIsoCheck(const float *x, int istart, int pplen, float ISO_LIMIT);

void pWaveRoutine(float *Xpc, int *Ponset, int *Poffset, int *Ppos);

int threeDataMax(float pWeight1, float pWeight2, float pWeight3);

void ECGMaxPosition(const float *sourceBeat, int length, float *returnMax, int *position);

void bufferCopy(const float *sourceBeat, int start, int finish, float *returnBuffer);

int localMaxAlgorithm(float *array, int begin, int finish);

void pSumCalculate(float *ecgdata, int Ponset, int Poffset, float *Psum);

void localMaxAlgorithm2(float *array, int begin, int finish, float max, int maxPosition, float min, int minPosition);

void LowPassFilterP(const float *X, float *Xpc);

void PWaveletDetect(const float *beat, unsigned int onset, int rr, float isoelectricLevel, int *Ponset, float *pPeak,
                    int *Poffset, int *pPeakPositionFirst, int *pPeakPositionSecond, int *type, float *Psum, int init);






//private:

//pwav_data_t2 result;
//pwav_data_t Presult[6];
//float  swa[5][500];
//float  swd[5][500];



//};




















#endif // PWAVEDETECTED_H
