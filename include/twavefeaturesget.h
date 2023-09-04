#ifndef TWAVEFEATURESGET_H
#define TWAVEFEATURESGET_H


typedef struct {
    int Tonset;
    int tOffset;
    int qtInterval;
    int type;
    float tAera;
    int tFirstPeak;
    int tSecondPeak;
    float tHeight;
} twav_data_t;


//class TwaveFeaturesGet
//{
//public:
//    TwaveFeaturesGet();
//    ~TwaveFeaturesGet();

//public:

int isoCheckT(const float *x, int istart, int len, float ISO_LIMIT);

//    [ max,maxPosition,min,minPosition ] = localMaxAlgorithm( array,begin,finish )
void
localMaxAlgorithmT(float *array, int begin, int finish, float *tmax, int *maxPosition, float *tmin, int *minPosition);

int IsoCheckDiffT(float *data, int j, int isoLength, float tlimit);

int GetMinimmumAngleT(const float *beat, int leftmost, int rightmost);

void LowPassFilterT(const float *X, float *Xpc);
//  void  tWavesFeaturesGet(float*beat,int onset,int rr,int offset,int SPeak,float isoelectricLevel ,float tLimit,float Crit,int*tFirstPeak,int*tSecondPeak,int*type,float*tHeight,int*Tonset,int*tOffset,float*tSum);

void tWavesFeaturesGet(float *beat, int onset, int rr, int offset, int SPeak, float isoelectricLevel, float tLimit,
                       float Crit,
                       twav_data_t *tPoint);//,int*tFirstPeak,int*tSecondPeak,int*type,float*tHeight,int*Tonset,int*tOffset,float*tSum);



//};

#endif // TWAVEFEATURESGET_H
