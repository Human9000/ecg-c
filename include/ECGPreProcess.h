#ifndef ECGPREPROCESS_H
#define ECGPREPROCESS_H


#define sampleRate          500
#define multiLeadNumberMax  12
#define cycleLength         10
#define iirFilterLength     7
#define lenghtLevelOne      51
#define lengthLevelTwo      351
#define constValue          403
#define ChannelMax          12


typedef  struct {
    float IIRpy[iirFilterLength];
    float IIRpx[iirFilterLength];
}IIRBandFilter_CTX;

typedef struct {
    float resultOne[lenghtLevelOne];
    float templateOne[lenghtLevelOne];
    float resultTwo[lengthLevelTwo];
    float templateTwo[lengthLevelTwo];
    float realTimeData[(lengthLevelTwo - 1) / 2];
    int count;
} baseLineRemoval_CTX;

typedef struct {
    float resultOneMulti[multiLeadNumberMax][lenghtLevelOne];
    float templateOneMulti[multiLeadNumberMax][lenghtLevelOne];
    float resultTwoMulti[multiLeadNumberMax][lengthLevelTwo];
    float templateTwoMulti[multiLeadNumberMax][lengthLevelTwo];
    float realTimeDataMulti[multiLeadNumberMax][(lengthLevelTwo - 1) / 2];
    int countMulti[12];
    float resultOne[lenghtLevelOne];
    float templateOne[lenghtLevelOne];
    float resultTwo[lengthLevelTwo];
    float templateTwo[lengthLevelTwo];
    float realTimeData[(lengthLevelTwo - 1) / 2];
    int count;
} baseLineRemovalMulti_CTX;

typedef struct {
    int flag;
    int c_threshold;
    float notchQueue[ChannelMax][cycleLength + 2];
    float notchData[ChannelMax][cycleLength];
    float pbuffer[ChannelMax][cycleLength];
    int index[ChannelMax];
} NotchFilter_50hz_CTX;
// float IntegerLowPass(float inputdata);

float baseLineRemoval(baseLineRemoval_CTX * ctx,  float inputdata, int initial) ;

void binaryInsertSort(float *resultOne, float currentData, float data, int lenght);

void dataSort(const float *inputArray, float *outputArray, int length);

float IIRBandFilter(IIRBandFilter_CTX *ctx,float inputdata, int initial);

float NotchFilter_50hz(NotchFilter_50hz_CTX *ctx,float datum, int channelNumber, int init);

float baseLineRemovalMulti(baseLineRemovalMulti_CTX *ctx,float ecgSample, int channel, int init);


#endif