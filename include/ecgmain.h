#ifndef ECGMAIN_H
#define ECGMAIN_H


#include "common_util.h"
#include "utility.h"
#include "ECGPreProcess.h"
#include "beatFeatureGet.h"
#include "twavefeaturesget.h"
#include "pWaveDetected.h"
#include "qrsDetect20230211.h"


#define   vTachyLength       3
#define   BEAT_MS90          45
//#define   sampleRate         500
#define   BEAT_MS1500        750
#define   BEAT_MS1600        800
#define   ECG_BUFFER_LENGTH  3000//1500
//#define   qrsFeaturesLength  500
#define   histLength         64
#define   templateLength     4
#define   ectopicArrayLength 7
#define   maxInteger         32767
#define   FIDMARK            200
#define   BEAT_MS200         100
#define   BEAT_MS220         110
#define   BEAT_MS140         70
#define   BEAT_MS130         65
#define   BEAT_MS50          25
#define   BEAT_MS800         400
#define   MS60               30
#define   BEAT_MS110         55
#define   BEAT_MS1000        500
#define   BEAT_MS500         250
#define   BEAT_MS3000        1500
#define   BEAT_MS4000        2000
#define   BEAT_SAMPLE_RATE   sampleRate
#define   ChannelMax          12
#define   vfArrayLength       6
//
#define  CoefficentConst   1//80//40//20//200//20//80//200


typedef enum {
    T_NORMAL = 0,                //normal
    T_INVERTED = 1,                //invert
    T_UPWARDS = 2,                //decreases monotonically
    T_DOWNWARDS = 3,            //increase  monotonically
} T_SHAPE;

//typedef enum  pType_main{
//   type_0,   type_1,   type_2,type_3,type_4,type_5,type_6,
//}P_SHAPE_main;

typedef struct ecgMainwdata_t {
    int rr;
    int onset;
    int offset;
    int qrsWidth;
    int Ponset;
    int Poffset;
    int Pwidth;
    float qrsHeight;
    float pPeak;
    int pType;
    int stChange;
    int stType;
    float speak;
    int noiseFlag;
    float tHeight;
    int Tonset;
    int tOffset;
    int qtInterval;
    float area;
    int prInterval;
    int ppInterval;
    float qrsMainWave;
    float tAera;
    float Psum;
    int tFirstPeak;
    int tSecondPeak;
    int rrDiff;
    int prDiff;
    int rrMax;
    float shannon;
    float ration;
    int noiseCount;
    int belowMs200;
    int morphologyType;
    float RRI;
    int type;
    int pPeakPositionFirst;
    int pPeakPositionSecond;
    float qPeak;
    int notchFlag;
    int vfFlag;

    int qpeak;
    //int     notchFlag;

} ecgMainwdata_t; // qrsFeatures

//supraventricular
typedef struct { // 房性早搏（APC）
    int apcTotal;                 // apc总计；
    int japcTotal;                // japc总计；
    int doubleApc;                // 双Apc；
    int apcBigeminyCount;         // apc二元计数；
    int apcTrigeminyCount;        // apc三叉计数；
    int atrialTachycardiaConut;   // 心房性心动过速；
    int longestApcTachycardia;    // 最长的Apc心动过速；
    int fastestApcTachy;          // 最快的Apc心动过速；
    int atrial_fibrillation;      // 心房颤动；
    int atrial_flutter;           // 心房扑动；
} supraventricular_type; // 室上型

//ventricular rhythm
typedef struct { // apc 室性早搏
    int pvcCount;                         // Pvc计数；
    int doublePvc;                        // 双Pvc；
    int vpcBigeminyCount;                 // vpcBigeminy计数；
    int vpcTrigeminyCount;                // vpc三叉神经计数；
    int ventricularTachycardiaCount;      // 室性心动过速计数；
    int longestPvcTachycardia;            // 最长Pvc心动过速；
    int fastestPvcTachy;                  // 最快的Pvc塔奇；
    int lbbbCount;                        // lbbb计数；
    int rbbbCount;                        // rbbb计数；

} ventricular_type; // 心室型 * 9

//sinus rhythm
typedef struct {
    int sinusArrhythmiaType;           // 窦性心律失常型；
    int sinusIrregularityType;         // 窦不规则型；
    int sinusArrestType;               // 窦性阻滞型；
    int sinusBradycardiaType;          // 窦性心动过缓型；
    int sinusTachycardiaType;          // 窦性心动过速型；
    int sickSinusSyndromeType;         // 病态窦房结综合征类型；

} sinus_type; // 窦型 * 6

//
typedef struct {
    int junctional_escape;
    int ventricular_escape;

} escape_type; // 逃逸类型 *2

//
typedef struct {

    int one_degreeCount;
    int two_degreeCount;
    int three_degreeCount;
    int high_degreeCount;
    int AsystolyCount;
} blocked_type; // 阻塞型 *5

//
typedef struct {

    float st_value;
    int type;

} st_segment;

//
typedef struct {
    supraventricular_type supra_arrhythmia; // 室上型 * 10
    ventricular_type ventricular_arrhythmia; // 心室型 * 9
    sinus_type sinus_arrhythmia; // 窦型 * 6
    escape_type escape_Count; // 逃逸类型 * 2
    blocked_type block_Count; // 阻塞型 * 5
    st_segment st_value; // ST段 * 2
} arrhythmia_type; // 心律失常型 * 34

void init_arrhythmia_type(arrhythmia_type *at);


// QRS波结构体
typedef struct {
    long q_point;            /*q点的相对于数据起点位置	单位:point  不存在则为-1*/
    long r_point;            /*r点的相对与数据起点位置	单位:point  不存在则为-1*/
    long j_point;            /*j点的相对于数据起点位置	单位:point  不存在则为-1*/
    int r_value;             /*r波的峰值 单位:unit*/
    int r_polar;             /*r波极性，TRUE波形朝上，FALSE波形朝下，ERROR未检测到极性*/
    char r_code;             /*qrs波类型*/
} ECG_QRSWave;

//T波结构体
typedef struct {
    int tw_value;            /* T波的峰值 单位:unit*/
    int t_polar;             /* T波极性，1波形朝上，0波形朝下，-1未检测到极性*/
    long ton_point;          /* T波起点的相对于数据起点位置	单位:point 不存在则为-1*/
    long toff_point;         /* T波终点的相对于数据起点位置	单位:point 不存在则为-1*/
    long tw_point;           /* T波峰点的相对于数据起点位置	单位:point 不存在则为-1*/
    T_SHAPE t_shape;

} ECG_TWave;


//心搏结构体
typedef struct {
    int id;                 /*标记心搏顺序号，从1开始*/
    int st_value[12];       /*st段的值  单位:unit*/
    int rr_interval;        /*rr间期 单位:unit 不存在则为-1*/
    int qrs_interval;       /*qrs时限，单位:unit 不存在则为-1*/
    int qrs_area;           /*qrs面积，单位:unit * unit 不存在则为-1*/
    int qt_interval;        /*qt间期 单位:unit 不存在则为-1*/
    int st_segment;         /*st段 单位:unit  不存在则为-1*/
    int hr;                 /*瞬时心率，单位：bpm*/
    int type;               //无意义
    char r_code;            /*r波码值，N代表正常，S代表室上性，V代表室性，X代表伪差，?未分类（疑问）不存在为标记为'*'*/
    ECG_QRSWave qrswave;    /*QRS波结构体*/
    ECG_TWave twave;        /*T波结构体*/
} ECG_BeatInfo;


//伪差结构体
typedef struct {
    int begin_point;      /*开始时间  单位:unit*/
    int end_point;        /*结束时间  单位:unit*/
} ECG_XInfo;


//心电数据结构体
typedef struct {
    unsigned char *tag;       //心电数据标记
    unsigned short *leadoff;  //导联脱落
    int *i;                   //I导联心电数据
    int *ii;                  //II导联心电数据
    int *iii;                 //III导联心电数据
    int *avR;                 //avR导联心电数据
    int *avL;                 //avL导联心电数据
    int *avF;                 //avF导联心电数据
    int *v;                   //V导联心电数据
} EcgData;


//HRV频域结构体
typedef struct {
    double ulf;                                   //超低频
    double vlf;                                   //极低频
    double lf;                                    //低频
    double hf;                                    //高频
    double lf_hf;                                 //低频与高频成分之比
    double tp;                                    //计算总功
} SimpleHRVFrequencyDomain;


//HRV时域结构体
typedef struct {
    float sdnn;                                //SDNN
    float sdann;                            //SDANN
    float sdannindex;                        //SDANNIndex
    float rmssd;                            //RMSSD
    float pnn50;                            //PNN50
    float msd;                                //MSD
    int nn50;                                //NN50
    float sdsd;                                //SDSD
} SimpleHRVTimeDomain;


typedef struct {
    int ri;                    //第一个标记节点
    int len;                //长度
} ECG_ARRNODE;

//房扑/房颤结构体
typedef struct {
    int begin;                //起始位置
    int end;                //结束为止
} ECG_AFIB;

//心律失常结构体
typedef struct {
    int bradynum;                                //心动过缓个数
    int *bradyarr;                                //心动过缓位置数组 外部malloc 数组
    int tachynum;                                //心动过速个数
    int *tachyarr;                                //心动过速位置数组 外部malloc 数组
    int pausenum;                                //长RR（停搏）个数
    int *pausearr;                                //长RR（停搏）位置数组 外部malloc 数组

    int vesum;                                    //室性早搏总数
    int *vesumarr;                                //室性早搏位置数组 外部malloc 数组
    int verun;                                    //室速个数
    ECG_ARRNODE *verunarr;                        //室速信息数组 外部malloc 数组
    int vetrigeminy;                            //室性三联律个数
    ECG_ARRNODE *vetrigeminyarr;                //室性三联律信息数组 外部malloc 数组
    int vebigeminy;                                //室性二联律个数
    ECG_ARRNODE *vebigeminyarr;                    //室性二联律信息数组 外部malloc 数组
    int vecouplet;                                //室性成对个数
    int *vecoupletarr;                            //室性成对信息数组 外部malloc 数组
    int vesingle;                                //室早单发个数
    int *vesinglearr;                            //室早单发位置数组 外部malloc 数组
    int velongestrun;                            //最长室速阵数
    int velongestrunnum;                        //最长室速位置

    int svesum;                                    //室上性早搏总数
    int *svesumarr;                                //室上性早搏位置数组 外部malloc 数组
    int sverun;                                    //室上速个数
    ECG_ARRNODE *sverunarr;                        //室上速信息数组 外部malloc 数组
    int svetrigeminy;                            //室上性三联律个数
    ECG_ARRNODE *svetrigeminyarr;                //室上性三联律信息数组 外部malloc 数组
    int svebigeminy;                            //室上性二联律个数
    ECG_ARRNODE *svebigeminyarr;                //室上性二联律信息数组 外部malloc 数组
    int svecouplet;                                //室上性成对个数
    int *svecoupletarr;                            //室上性成对信息数组 外部malloc 数组
    int svesingle;                                //室上早单发个数
    int *svesinglearr;                            //室上早单发位置数组 外部malloc 数组

    int svelongestrun;                            //最长室上速阵数
    int svelongestrunnum;                        //最长室上速位置


    int threeAB;                                //房室三度传导阻滞个数
    int *threeABarr;                            //房室三度传导阻滞位置   外部malloc 数组
    int towAB;                                    //房室二度传导阻滞个数
    int *towABarr;                                //房室二度传导阻滞位置   外部malloc 数组

    ECG_AFIB *af;                                //房扑/房颤位置信息
    int afnum;                                    //房扑/房颤个数
} ECG_ARRHYTHMIA;


typedef struct {
    int queueLength;
    float Sdeep[10];
    float informArray[8][100];
} PvcLibInfor;

typedef struct {
    float qrsArea;
    int qrsWidth;
} PvcNodeStru;

// static
typedef struct {
    int sinusIrregularityConut;
    int sinusArrestCount;
    int sinusBradycardiaCount;
    int sinusTachycardiaCount;
    int sickSinusSyndromeCount;
    int modeCycle;
    int sinoAtrialArray[10];
} sinusArrhythmia_CTX;

typedef struct {
    int AVArray[10];
    int firstDegreeBlockCount;
} firstDegreeAVBlock_CTX;

typedef struct {
    int vpcString[vTachyLength];
    int vpcBigeminy;
    int vpcTrigeminy;
    int vpcBigString[vTachyLength + 1];
    int vpcTriString[vTachyLength * 2];
    int longestPvcTachycardiaCount;
    int pvcTachycardiaSegment;
    int pvcPointerCount;
    int pvcBiSegmentCount;
    int pvcBiLongest;
    int fastestPvcTachycardia;
    int doublePvcCount;
    int doubleArray[vTachyLength];
    int rrSum;
    int rrString[vTachyLength];
    int biFlag;
    int flag;
    int count;
    int triCount;
    int triFisrFlag;
    int triSecongFlag;
    int pvcTriSegmentCount;
    int pvcTriLongest;
} ventricularPrematureCount_CTX;


typedef struct {
    int apcString[vTachyLength];
    int apcBigeminy;
    int apcTrigeminy;
    int apcBigString[vTachyLength + 1];
    int apcTriString[vTachyLength * 2];
    int longestApcTachycardiaCount;
    int apcTachycardiaSegment;
    int apcPointerCount;
    int apcBiSegmentCount;
    int apcBiLongest;
    int fastestApcTachycardia;
    int doubleApcCount;
    int apcdoubleArray[vTachyLength];
    int apcrrSum;
    int apcrrString[vTachyLength];
    int apcbiFlag;
    int apcflag;
    int apcCount;
    int apctriCount;
    int apctriFisrFlag;
    int apctriSecongFlag;
    int apcTriSegmentCount;
    int apcTriLongest;
    float normalTemplate[BEAT_MS200];
} atrialPrematureCount_CTX;


typedef struct {
    int noiseFlagArray[3];
    int noiseTag;
} QRSMorphology_CTX;

typedef struct {
    BeatFeaturesGet_CTX beatFeaturesGetCtx;
} multiLeadsFeaturesGet_CTX;

typedef struct {
    float tempBeat[sampleRate];
    float newArray[BEAT_MS200];
    int pvcNumber;
    float ectopicArray[ectopicArrayLength];
    float tempectopicArray[ectopicArrayLength];
    float rrAFArray[histLength];
    int ECGBufferIndex;
    float ECGBuffer[ChannelMax][ECG_BUFFER_LENGTH];
    float BufferV1[ECG_BUFFER_LENGTH];
    float rWaveLocation[ECG_BUFFER_LENGTH / CoefficentConst];
    float rWaveLocationV1[ECG_BUFFER_LENGTH / CoefficentConst];
    float lastBeat[sampleRate];
    float secondBeat[sampleRate];
    float prDistribution[histLength];
    int lastRR;
    int lastRRV1;
    int peakNumber;
    int indexFlag;
    int learningNumber;
    int learningFlag;
    int complementFlag;
    int afEttopicCount;
    int rrAFCount;
    int pvcAFlag;
    int apcAFlag;
    int apcAFlagCount;
    int belowMs200;
    int avrQRSRR;
    int noiseCount;
    int ppIntervalTotal;
    int ppIntervalCount;
    float currentAverSpeak;
    float maxQRSAerea;
    float ECGBufferMulti[ChannelMax][ECG_BUFFER_LENGTH / CoefficentConst];
    int globalIndex[ChannelMax];
    int rbbbCount;
    int lbbbCount;
    int japc;
    int apc;
    int wpwCount;
    int doublePvc;
    int vpcBigeminyCount;
    int vpcTrigeminyCount;
    int ventricularTachycardiaCount;
    int longestPvcTachycardia;
    int fastestPvcTachy;
    int doubleAvc;
    int apcBigeminyCount;
    int atrialTachycardiaConut;
    int longestApcTachycardia;
    int apcTrigeminyCount;
    int fastestApcTachy;
    int ventrEscapeCount;
    int junctEscapeCount;
    int templateFlag;
    int sinusArrhythmiaType;
    int sinusIrregularityType;
    int sinusArrestType;
    int sinusTachycardiaType;
    int sickSinusSyndromeType;
    int sinusBradycardiaType;
    int twoDgreeBlockNumber;
    float avrQRSWidth;
    float currentTArea;
    float currentArea;
    float currentPeak;
    float currentPsum;
    float tHeightAverage;
    float currentRRSMD;
    float qrsRR[templateLength];
    float qrsWidth[templateLength];
    float qrsSpeak[templateLength];
    float qrsArea[templateLength];
    float tWaveArea[templateLength];
    float Pheight[templateLength];
    float Parea[templateLength];
    float tNormalAverage[templateLength];
    float pwaveTemplate[2 * BEAT_MS50];
    float apcProcessing[MS60];
    float AFapcProcessing[MS60];
    PvcLibInfor pvcLibInfo;
    int stChangeCount;
    int pvcLength;
    int currentWidthMax;
    int currentQRS;
    float pvcSum;
    float pvcStandard;
    float pvcWidthAverage;
    int pvc;
    PvcNodeStru pvcArray[templateLength];
    float pvcStandardDerivation;
    float pvcStandardWidth;
    float pvcWidthMaximum;
    float pvcAreaMaximum;
    int pvcTemplateCount;
    int ST_IschaemicEpisode[20];
    int pvcWidthArray[6];
    int pvcQrsWidthArray[6];
    int effecitiveAFcount;
    int qrsFeaturesCountPreSec;
    int qrsFeaturesCountPre;
    int qrsFeaturesCount;
    ecgMainwdata_t qrsFeatures[ECG_BUFFER_LENGTH / CoefficentConst];
    float normalTemplate[BEAT_MS200];
    int afEntryFlag;
    int reEntryFlag_II;
    unsigned int reEctopicFlag;
    int lastValue;
    int rrCountP;
    int indexRr;
    int rrInteral[ECG_BUFFER_LENGTH];
    ecgMainwdata_t MultiQrstFeatures[ECG_BUFFER_LENGTH];
    int MultiPeakNumber;
    int n;
    int sIndex;
    int rrWaveLocation[ECG_BUFFER_LENGTH / CoefficentConst];
    int testCount;
    IntegerLowPassMulti_CTX integerLowPassMultiCtx;
    sinusArrhythmia_CTX sinusArrhythmiaCtx;
    firstDegreeAVBlock_CTX firstDegreeAvBlockCtx;
    ventricularPrematureCount_CTX ventricularPrematureCountCtx;
    atrialPrematureCount_CTX atrialPrematureCountCtx;
    QRSMorphology_CTX qrsMorphologyCtx;
    IIRBandFilter_CTX iirBandFilterCtx;
    baseLineRemoval_CTX baseLineRemovalCtx;
    IntegerLowPass_CTX integerLowPassCtx;
    qrsDetectMultiLeads_CTX qrsDetectMultiLeadsCtx;
    baseLineRemovalMulti_CTX baseLineRemovalMultiCtx;
    NotchFilter_50hz_CTX notchFilter50HzCtx;
    multiLeadsFeaturesGet_CTX multiLeadsFeaturesGetCtx;


    int count;
    int rrCount;
    arrhythmia_type ret;

} ArrhythmiaAnalysis20230311_CTX;
//class ecgMain
//{
//public:
//    ecgMain();
//    ~ ecgMain();
//public:



int wpwCalculate(float *beat, int onset, float qrsHeight);

float FeaturesOptimize(float currentFeature, float currentValue);

void templateUpdate(PvcLibInfor *array, float *newArray);

float arrayCreate(float *lastBeat, float *tempBeat, int posiFirst, int posSecond, int len);

float calculateCrit(const float *beat, int rrValue);

int findSinusPWaves(float *array, int rr, int Interval, int currentOnset, int tOffset, float thresh, float currentPsum);

void findKnots(float *y, int *knotsPosition, int *originalPosition, float *firstDiffer);


void sinusArrhythmia(sinusArrhythmia_CTX *CTX, int ppInterval, int *sinusArrhythmiaType, int *sinusIrregularityType,
                     int *sinusArrestType,
                     int *sinusBradycardiaType, int *sinusTachycardiaType, int *sickSinusSyndromeType, int init);

//void sinusArrhythmia(int ppInterval, int *sinusArrhythmiaType, int *sinusIrregularityType, int *sinusArrestType,
//                     int *sinusBradycardiaType, int *sinusTachycardiaType, int *sickSinusSyndromeType, int init);

int firstDegreeAVBlock(firstDegreeAVBlock_CTX *ctx, int prInterval, int init);

void
ST_Morphology(const float *array, int begin, int finish, float amplifierGainConst, int len, int tFirstPeak, float Crit,
              int *type);

void ST_SegmentAnalysis(float *array, int rr, int pStart, int pOffset, int onset, int offset, float amplifierGainConst,
                        int tFirstPeak, float Crit, float leadStandard, int *ST_MorphologyType, int *ST_Change,
                        float *ST_ElavatorAmph, int *type);

void ventricularDataCopy(int doublePvcCount, int pvcBiSegmentCount, int pvcTriSegmentCount, int pvcTachycardiaSegment,
                         int longestPvcTachycardiaCount, int fastestPvcTachycardia, int *doublePvc,
                         int *vpcBigeminyCount, int *vpcTrigeminyCount, int *ventricularTachycardiaConut,
                         int *longestPvcTachycardia, int *fastestPvcTachy);

void clearBuffer(int *vpcTriString, int *doubleArray, int *vpcBigString);

void
ventricularPrematureCount(ventricularPrematureCount_CTX *ctx, int pvcFlag, int rrInterval, int *doublePvc,
                          int *vpcBigeminyCount, int *vpcTrigeminyCount,
                          int *ventricularTachycardiaCount, int *longestPvcTachycardia, int *fastestPvcTachy, int init);

float templateSeek(PvcLibInfor *array, float *source);

float shannonDistribution(float *rrArray, int len);

float correlationCoefficient(const float *array1, const float *array2, int len);

float FuzzyReasoning(int classification, float inputData);

float templateCompute(const float *array, int number);

void arrayCopy(const float *sourceArray, float *destArray, int len);

void findTruePeak(float *array, int rrValue, float *adjustedPeak);

//void QRSMorphology(float *beat, int flag, int *noiseFlag, int *morphologyType, int init);
void QRSMorphology(QRSMorphology_CTX *ctx, float *beat, int flag, int *noiseFlag, int *morphologyType, int init);

int rWavePeak(float *beat, int rrValue, float Crit);

//    int   counterAdjust(int counter, int modeValue );  // 此函数已经移动到 common_util.h 中
void DataSort(float *dataArray, int length);

void
shannonEntropyCompute(float *rrArray, const float *prDistribution, float *shannonEntropy, float *RMSSD,
                      float *prIntervalDiff,
                      float *RRImaxCount);

void atrialPrematureCount(atrialPrematureCount_CTX *ctx, int apcFlag, int rr, int *doubleAvc, int *apcBigeminyCount,
                          int *apcTrigeminyCount,
                          int *atrialTachycardiaConut, int *longestApcTachycardia, int *fastestApcTachy, int init);

int ArrhythmiaAnalysis(float IILeadData, float V1LeadData, float V5LeadData, int *rrCount, arrhythmia_type *result,
                       int init);

void templateCopy(const float *sourceBeat, float *destData);

//int ArrhythmiaAnalysis20230311(const float *array, int ChannelNo, int *rrCount, arrhythmia_type *result, int init,
//                               int test_data_no);

int ArrhythmiaAnalysis20230311(ArrhythmiaAnalysis20230311_CTX *ctx, const float *array, int ChannelNo, int *rrCount,
                               arrhythmia_type *result, int init,
                               int test_data_no);

void
multiLeadsFeaturesGet(multiLeadsFeaturesGet_CTX *ctx, float ECGBuffer[][ECG_BUFFER_LENGTH], int rrValue, int *onset,
                      int *offset, int *QPeak,
                      float *speakheight, float *rweight, int *beatBegin, int *beatEnd, float *amp, int *noiseFlag,
                      float *isoelectricLevel, int *notchFlag, float *IIArea, int temp);

//private:

//  static  PvcLibInfor pvcLibInfo;
//  float  newArray[BEAT_MS200];
//  float  ectopicArray[ectopicArrayLength];
//  float  tempectopicArray[ectopicArrayLength];
//  //int    nullpointer[4];
//  int    pvcNumber;
//  float  morphologyIndex;
//  float  recentPrInterval;
//  //
//  int    ST_MorphologyType;
//  int    ST_Change;
//  float  ST_ElavatorAmph;
//  int    type;
//  float  amplifierGainConst;
//  float  leadStandard;
//  float  fuzzyNormalResult;
//  float  areaIndex;
//  float  fuzzyPVCResult;
//  float  fuzzyApcResult;

////
//  static int stChangeCount;
//  static int ST_IschaemicEpisode[20];
//  static int noiseFlagArray[3];
//  static int noiseTag;

// //
//  static int  sinusArrhythmiaType;
//  static int  sinusIrregularityType;
//  static int  sinusArrestType;
//  static int  sinusTachycardiaType;
//  static int  sickSinusSyndromeType;
//  static int  sinusBradycardiaType;

//  //
//  static int sinusIrregularityConut;
//  static int sinusArrestCount;
//  static int sinusBradycardiaCount;
//  static int sinusTachycardiaCount;
//  static int sickSinusSyndromeCount;

//  static int modeCycle;
//  static int templateFlag;
//  static int sinoAtrialArray[10];

////for first degree block
//static int AVArray[10];
//static int firstDegreeBlockCount;

//// for ventricularPrematureCount
// static int  vpcString[vTachyLength];
// static int  vpcBigeminy;
// static int  vpcTrigeminy;
// static int  vpcBigString[vTachyLength+1];
// static int  vpcTriString[vTachyLength*2];
// static int  longestPvcTachycardiaCount;
// static int  pvcTachycardiaSegment;
// static int  pvcPointerCount;
// static int  pvcBiSegmentCount;
// static int  pvcBiLongest;
// static int  fastestPvcTachycardia;
// static int  doublePvcCount;
// static int  doubleArray[vTachyLength];
// static int  rrSum;
// static int  rrString[vTachyLength];
// static int  biFlag;
// static int  flag;
// static int  count;
// static int  triCount;
// static int  triFisrFlag;
// static int  triSecongFlag;
// static int  pvcTriSegmentCount;
// static int  pvcTriLongest;

// static float pvcStandardDerivation;
// static float pvcStandardWidth;
// static float pvcWidthMaximum;
// static float pvcAreaMaximum;
// static int   pvcTemplateCount;

//  static int   testCount;/////////////////////////////////////////////////
// //for void ecgMain:: atrialPrematureCount
// static int  apcString[vTachyLength];
// static int  apcBigeminy;
// static int  apcTrigeminy;
// static int  apcBigString[vTachyLength+1];
// static int  apcTriString[vTachyLength*2];
// static int  longestApcTachycardiaCount;
// static int  apcTachycardiaSegment;
// static int  apcPointerCount;
// static int  apcBiSegmentCount;
// static int  apcBiLongest;
// static int  fastestApcTachycardia;
// static int  doubleApcCount;
// static int  apcdoubleArray[vTachyLength];
// static int  apcrrSum;
// static int  apcrrString[vTachyLength];
// static int  apcbiFlag;
// static int  apcflag;
// static int  apcCount;
// static int  apctriCount;
// static int  apctriFisrFlag;
// static int  apctriSecongFlag;
// static int  apcTriSegmentCount;
// static int  apcTriLongest;

////
// static  int japc;
// static  int apc;
// static  int pvc;
// static  int wpwCount;
// static int rbbbCount;
// static int lbbbCount;

// //
// static int doublePvc;
// static int vpcBigeminyCount;
// static int vpcTrigeminyCount;
// static int ventricularTachycardiaCount;
// static int longestPvcTachycardia;
// static int fastestPvcTachy;

//  //
// static int doubleAvc;
// static int apcBigeminyCount;
// static int atrialTachycardiaConut;
// static int longestApcTachycardia;
// static int apcTrigeminyCount;
// static int fastestApcTachy;

// static int ventrEscapeCount;
// static int junctEscapeCount;
// static int twoDgreeBlockNumber;


// static int ECGBufferIndex;
// static float ECGBuffer[ECG_BUFFER_LENGTH];
// static float BufferV1[ECG_BUFFER_LENGTH];


// static int lastRR;
// static int lastRRV1;
// static int peakNumber;
// static int peakNumberV1;
// static int indexFlag;

// static float lastBeat[sampleRate];
// static float secondBeat[sampleRate];
// static float rWaveLocation[ECG_BUFFER_LENGTH];
// static float rWaveLocationV1[ECG_BUFFER_LENGTH];
// static float  rrAFArray[histLength];
// static  float  prDistribution[histLength];


// static int learningNumber;
// static int learningFlag;
// static int complementFlag;
// static int afEttopicCount;
// static int rrAFCount;
// static int pvcAFlag;
// static int apcAFlag;
// static int apcAFlagCount;
// static int belowMs200;
// static int avrQRSRR;
// static int noiseCount;
// static float currentAverSpeak;
// static int pvcLength;
// static int currentWidthMax;
// static int currentQRS;

// static int  ppIntervalTotal;
// static int  ppIntervalCount;

// // static float  avrQRSRR;
// static float  avrQRSWidth;
// static float  currentTArea;
// static float  currentArea;
// static float  currentPeak;
// static float  currentPsum;
// static float  tHeightAverage;
// static float  maxQRSAerea;

// static  float currentRRSMD;
// static int effecitiveAFcount;


// static  float pvcSum;
// static  float pvcStandard;
// static  float pvcWidthAverage;

// float   tempBeat[sampleRate];
// static ecgMainwdata_t qrsFeatures[ECG_BUFFER_LENGTH];

// static  float qrsRR[templateLength];//=zeros(1,4);
// static  float qrsWidth[templateLength];//=zeros(1,4);
// static  float qrsSpeak[templateLength];//=zeros(1,4);
// static  float qrsArea[templateLength];//=zeros(1,4);
// static  float tWaveArea[templateLength];//=zeros(1,4);
// static  float Pheight[templateLength];//=zeros(1,4);
// static  float Parea[templateLength];//=zeros(1,4);
// static  float tNormalAverage[templateLength];//=zeros(1,4);

// static int pvcWidthArray[6];
// static int pvcQrsWidthArray[6];

// static PvcNodeStru pvcArray[templateLength];

// static  float normalTemplate[BEAT_MS200];
// static  float  pwaveTemplate[2*BEAT_MS50];

// static  float apcProcessing[MS60 ];
// static  float AFapcProcessing[MS60 ];

//};

#endif // ECGMAIN_H
