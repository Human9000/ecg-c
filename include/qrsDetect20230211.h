#ifndef QRSDETECT20230211_H
#define QRSDETECT20230211_H

#define EPS 0.000001
//%%%%%%%%%%     500HZ
#define MS800 400
#define MS1000 500
#define MS1500 750
#define MS180 90
#define MS150 75
#define MS340 170
#define MS360 180
#define MS280 140

#define MS95 48
#define MS80 40
#define MS200 100
#define MS120 60
#define MS140 70
#define MS290 145
#define MS100 50
#define MS400 200
#define BEAT_MS10 5
#define maxQueueLength 3000//%3000;%8000;%%1000000;  3000;%%
#define ModNumber maxQueueLength
#define QRSDelayLength 148
#define WINDOW_WIDTH 40
#define ForgetFactor 0.05
#define WeightFactor 0.15
#define PRE_BLANK 100
#define FILTER_DELAY (21+MS200)
#define fs   500
#define maxNegative (-32768)
#define MIN_PEAK_AMP 10
#define TH 0.2125
//#define  CoefficentConst   1
#define miniPeak 10

#define leadsNumber 12

//class qrsRoute
//{

//public:
//    qrsRoute();
//    ~qrsRoute();

//public:


typedef struct {
    int RRcount;
    int indexMulti[leadsNumber];
    float tempPeakMulti[leadsNumber];
    float x_dervMulti[leadsNumber][10];
    float MovingWindowSumMulti[leadsNumber];
    int MovingWindowSumPtrMulti[leadsNumber];
    float dataMulti[leadsNumber][40];
    int timeSinceMaxMulti[leadsNumber];
    float lastDatumMulti[leadsNumber];
    float maxMulti[leadsNumber];
    float MIN_PEAK_AMPMulti[leadsNumber];
    int preBlankCntMulti[leadsNumber];
    int qpkcntMulti[leadsNumber];
    int countMulti[leadsNumber];
    int initBlankMulti[leadsNumber];
    float rrbufMulti[leadsNumber][8];
    float noiseMulti[leadsNumber][8];
    float newPeakMulti[leadsNumber];
    float qrsbufMulti[leadsNumber][8];
    float rsetBuffMulti[leadsNumber][8];
    float pyMulti[leadsNumber][8];
    float pxMulti[leadsNumber][8];
    float det_threshMulti[leadsNumber];
    int rsetCountMulti[leadsNumber];
    float nmeanMulti[leadsNumber];
    float rrmeanMulti[leadsNumber];
    float qmeanMulti[leadsNumber];
    int sblocMulti[leadsNumber];
    int QrsDelayMulti[leadsNumber];
    float sbpeakMulti[leadsNumber];
    float initMaxMulti[leadsNumber];
    int lastPositionMulti[leadsNumber];
    int rMulti[leadsNumber];
    int qrscountMulti[leadsNumber][maxQueueLength];
    int adjcentPeakMulti[leadsNumber];
    float currentNewPeakMulti[leadsNumber];
    int lastNewPeakPositionMulti[leadsNumber];
    int currentNewPeakPositionMulti[leadsNumber];
    float ecgBufferMulti[leadsNumber][maxQueueLength];
    float totalBufferMulti[maxQueueLength];
    float PeakPositionMulti[leadsNumber][maxQueueLength];
    int r;
    float lostPeakTemp;
    int MultiVector[leadsNumber];
    int VectroPointer;
    int multiRRCount;
    int multiRRVector[maxQueueLength];
    int sbcount;
} qrsDetectMultiLeads_CTX;

int pqrsDetect2023(float inputdata, float *ECGBuffer, int init);//,int peakNumber,float*ECGBuffer);
float qrsHFNoiseCheck(float *beat, float qrsMin, float qrsMax);

int qrsDetectMultiLeads(qrsDetectMultiLeads_CTX *ctx,float inputdata, int multiLeadNumber, int init);
//int qrsDetect2023(float inputdata,int multiLeadNumber,float *ECGBuffer,int init);
//   int filtfilt( double* x, double* y, int xlen, double* a, double* b, int nfilt);
//  void filter(const double* x, double* y, int xlen, double* a, double* b, int nfilt);
//public:
//  // static float TH;
//   static int   index;
//   static float tempPeak ;
//   static float x_derv[10];
//   static float MovingWindowSum;
//   static int   MovingWindowSumPtr;
//   static float data[WINDOW_WIDTH];
//   static int   timeSinceMax;
//   static float lastDatum;
//   static float max;
// //  static float MIN_PEAK_AMP;
//   static int   preBlankCnt;
//   static int   qpkcnt;
//   static int   count;
//   static int initBlank;
//   static float rrbuf[8];
//   static float noise[8];
//   static float newPeak;
//   static float qrsbuf[8];
//   static float rsetBuff[8];
//   static float det_thresh;
//   static int   rsetCount;
//   static float nmean;
//   static float rrmean;
//   static float qmean;
//   static int   sbloc;
//   static int   QrsDelay;
//   static float sbpeak ;
//   static float initMax;
//   static float lastPeak;
//   static int   lastPeakIndex;
//   static int   currentPosition1;
//   static int   lastPosition;
//   static float PeakPosition[maxQueueLength];
//   static int   r;
//   static int   qrscount[maxQueueLength];
//   static float py[8];
//   static float px[8];
//   static int  adjcentPeak;
// //  static int   rrBuffer[maxQueueLength];
//   static int   RRcount;
//   static float qrsFeatures;
//   static float firstHeight;
//   static float secondHeight;
//   static float lastNewPeak;
//   static float currentNewPeak;
//   static int   lastNewPeakPosition;
//   static int   currentNewPeakPosition;
//   static int   smallCount;
// //  static int   buffer[maxQueueLength];
//   static float rrArray[5];
//   static float lostPeakArray[6];
//   static int   lostPeakPointer;
//   static float lostPeakTemp;
//   static int   lostPeakOutputPointer;
//   static float monintorArray[6];
//   static int   monitoringPointer;
//   static int   RRInterval;
//   static float ECGBuffer[maxQueueLength];
//   static int   sbcount;
//  // static float result2[500];



//};



#endif // QRSDETECT20230211_H
