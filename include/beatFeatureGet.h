#ifndef  BEATFEATUREGET_H
#define  BEATFEATUREGET_H

#define  sampleRate       500
#define  arrayLimit       300


//class qrsFeaturesClass
//{

//public:

//    qrsFeaturesClass();
//    ~qrsFeaturesClass();

//  public:


typedef struct {
    float noiseFlagArray[3];
    int noiseTag;
} BeatFeaturesGet_CTX;


float pIsoLevelB(const float *array, int position);

void qrsBoundaryB(const float *beat, float Crit, float mainWaveHeight, int *QRSbegin, int *QRSfinish, int *QRSPeaksPositions,
                  float *QRSPeaks, int *qrsIndex);

int IsoCheckB(const float *data, int j, int isoLength);

void LowPassFilterB(const float *X, int Fs, int Fpa, float *Xpc);

int IsoCheckDiffB(float *data, int j, int isoLength, float limit);

int IsoCheckDiffRightB(float *data, int j, int isoLength, float limit);

int GetMinimmumAngleB(const float *beat, int leftmost, int rightmost);

void BeatFeaturesGet(BeatFeaturesGet_CTX *ctx, float *tempBeat, int *SPeak, int interval, int *onset, int *offset, int *QPeak, float *speakheight,
                     float *rweight, int *noiseFlag, float *isoelectricLevel, int *notchFlag, int init);
// void  [onset,offset,QPeak,speakheight,rweight,beatBegin,beatEnd,amp,noiseFlag ,isoelectricLevel,notchFlag,IIArea] = multiLeadsFeaturesGet(ECGBuffer,rrValue,temp)

//private:

//    float  firstDiff[arrayLimit];
//    float  secondDiff[arrayLimit];
//    float  fnew[arrayLimit];
//    float  data[arrayLimit];
//    float  result[arrayLimit];
//    float  y[arrayLimit+5];

//    static int noiseFlagArray[3];
//    static int noiseTag;

//    float QRSPeaks[10];//=zeros(1,BEAT_MS20);
//    int   QRSPeaksPositions[10];//=zeros(1,BEAT_MS20);
//    float Xpc[sampleRate];//=zeros(1,dataLength);


//};













































#endif // BEATFEATUREGET_H
