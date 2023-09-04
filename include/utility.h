#ifndef UTILITY_H
#define UTILITY_H


#define multiLeadNumberMax 12

typedef struct {
    float data[10];
    float y0;
    float y1;
    float y2;
} IntegerLowPass_CTX;


typedef struct {
    double data[10];
    double y0;
    double y1;
    double y2;
    double dataMulti[multiLeadNumberMax][10];
    double y0Multi[multiLeadNumberMax];
    double y1Multi[multiLeadNumberMax];
    double y2Multi[multiLeadNumberMax];
    int mpa;
    int count;
}IntegerLowPassMulti_CTX;

//public:
int IsoCheckDiffRight(float *data, int j, int isoLength, float limit);

void ecgCopy(const float *array, float *result, int start, int finish);

float ecgMean(const float *array, int len);

void ecgDiff(const float *array, float *result, int len);

float noiseLevel(const float *array, int len);

float HFnoiseCheck(float *array);

int searchPeakRoute(const float *ECGBuffer, int searchindex);

float ECGMax(const float *array, int length);

float ECGMin(const float *array, int length);

float ecgSum(float *array, int len);

void LowPassFilter(const float *X, float *Xpc);


//realtime
// float  IntegerLowPass(float datum);
float IntegerLowPass(IntegerLowPass_CTX *ctx,float datum, int init);
// float  IntegerLowPassV1(float datum);
float IntegerLowPassV1(float datum, int init);


// MultiLeads
//float IntegerLowPassMulti(float datum, int multiLeadNumber, int init);
float IntegerLowPassMulti(IntegerLowPassMulti_CTX *ctx, float datum, int multiLeadNumber, int init);

// multiLeads

////private:
////  lead II
//   static float data[10];
//   static float y0;
//   static float y1;
//   static float y2;

////Lead V1
//   static float data_V1[10];
//   static float y0_V1;
//   static float y1_V1;
//   static float y2_V1;

//};

#endif // UTILITY_H
