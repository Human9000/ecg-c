#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common_util.h"
#include "utility.h"
#include "ECGPreProcess.h"
#include "beatFeatureGet.h"
#include "twavefeaturesget.h"
#include "pWaveDetected.h"
#include "qrsDetect20230211.h"
#include "ecgmain.h"


#define  BEAT_MS10     5
#define  BEAT_MS30    15
#define  ampValue     200

#define  pvcQueLength 6

#define  nomalRation  1
#define  shannLen     64
#define  BEATLGTH     500
#define  templateLength 4
#define  apcProcessLength 20// BEAT_MS60
#define  SAMPLE_RATE   500

#define  pvcThrehold   0.580
#define  Fs           sampleRate
#define  BEAT_MS8    4
#define  BEAT_MS80   40
#define  BEAT_MS20   10
#define  BEAT_MS160  80
#define  BEAT_MS60   30


#define  BEAT_MS250  125
#define  BEAT_MS40   20
#define  arrayLimit  300
#define  BEAT_MS480  240
#define  BEAT_MS6    3
#define  sinusArrLen 10

//#define   MS80        40
#define   BEAT_MS10   5


#define   BEAT_MS1200  600
#define   BEAT_MS120   60
#define   histWidth    16

#define  BEAT_MS64     32
#define  BEAT_MS72     36

#define  BEAT_MS600    300
#define  BEAT_MS1100   550
#define  BEAT_MS540    270

#define  BEAT_MS2000   1000
#define BEAT_MS230     115

#define  BEAT_MS100    50
#define  BEAT_MS400    200
#define  BEAT_MS360    180
#define  peakLocationMod  ECG_BUFFER_LENGTH
//#define  amplifierGainConst    200
//#define  leadStandard          0.2

#define BEAT_MS150     75
#define  maxInteger    32767

//#pragma warning(disable:4244)

//ecgMain::ecgMain()
//{

//    //tempBeat[sampleRate]={0};
//    memset(tempBeat,0,sampleRate*sizeof(float));
//    memset(rrAFArray,0,histLength*sizeof(float));
//    memset(newArray,0,BEAT_MS200*sizeof(float));
//    memset(ectopicArray,0,ectopicArrayLength*sizeof(float));
//    memset(tempectopicArray,0,ectopicArrayLength*sizeof(float));

//    recentPrInterval=0;
//    morphologyIndex=0;
//    pvcNumber=0;

//    ST_MorphologyType=0;
//    ST_Change=0;
//    ST_ElavatorAmph=0;
//    type=0;

//    amplifierGainConst =200;
//    leadStandard =0.2;
//    fuzzyNormalResult=0;
//    fuzzyApcResult=0;
//    areaIndex=0;
//    fuzzyPVCResult=0;
//   // nullpointer[4]={0};

//}

//ecgMain::~ecgMain()
//{

//}

//ecgMainwdata_t ecgMain::qrsFeatures[ECG_BUFFER_LENGTH]={{0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0}};





// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int wpwCalculate(float *beat, int onset, float qrsHeight) {

//    int BEAT_MS15=8;
//    int ampliferFactor=200;
    int pwpResult = 0;
    // % Find the izoelectric point on Average ecg
    int Flat = 2 * BEAT_MS10;       // % 20 ms
    int From = FIDMARK - Flat;
    int To = From - BEAT_MS80;
    float peak_y = qrsHeight;
    if (To < 10)
        To = 10;
    //end
//    float Crit=0.02*peak_y;	//% make Crit proportional to the QRSampl;



    return pwpResult;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
float FeaturesOptimize(float currentFeature, float currentValue) {
    float coef = currentValue / currentFeature;
    int weight;
    if ((coef >= 0.9) && (coef <= 1.1))
        weight = 1;
    else if (((coef >= 0.7) && (coef <= 0.9)) || ((coef >= 1.1) && (coef <= 1.5)))
        weight = 3;
    else if (((coef >= 0.4) && (coef <= 0.7)) || ((coef >= 1.5) && (coef <= 1.8)))
        weight = 7;
    else
        weight = 8;

    double optimizationValue = (currentFeature + currentValue * (float) weight) / (float) (weight + 1);
    return (float) optimizationValue;

}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void templateUpdate(PvcLibInfor *array, float *newArray) {
    int queMaxLength = 8;
    int len = 100;
    int maxLen = array->queueLength;
    float temp[100] = {0};//zeros(1,100);
    int num = 0;
    float maxCoef = 0;
    int maxLevel = 0;
    float correlation;

    int i;
    while (num < maxLen) {
        i = 0;
        while (i < 100) {
            temp[i] = array->informArray[num][i];
            i = i + 1;
        }
        correlation = correlationCoefficient(temp, newArray, len);// correlationCoefficient( temp,newArray );
        if (correlation > maxCoef) {
            maxCoef = correlation;
            maxLevel = num;
        }
        num = num + 1;
    }//end
    if (maxCoef > 0.92) {
        i = 0;
        while (i < 100) {
            array->informArray[maxLevel][i] = (float) (temp[i] + newArray[i]) / 2;
            i = i + 1;
        }
    } else if ((maxCoef < 0.92) && (maxLen < queMaxLength)) {
        i = 0;
        while (i < 100) {
            array->informArray[maxLen][i] = newArray[i];
            i = i + 1;
        }
        array->queueLength = array->queueLength + 1;
        maxLen = maxLen + 1;
    }

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float arrayCreate(float *lastBeat, float *tempBeat, int posiFirst, int posSecond, int len) {

    float array1[50];//=zeros(1,50);
    float array2[50];//=zeros(1,50);
    //arraySlope=zeros(1,20);
    int errorLength = 10;
    // %  global numberIndex;

    // %  temp=abs(lastBeat(posiFirst));
    int k = posiFirst - errorLength;
    float temp = fabsf(lastBeat[k]);
    int i = k + 1;
    while (i < (posiFirst + errorLength)) {
        if (fabsf(lastBeat[i]) > temp) {
            temp = fabsf(lastBeat[i]);
            k = i;
        }
        i = i + 1;
    }
    if (k != posiFirst)
        posiFirst = k;
    //end

    int kk = posiFirst - len;
    int mm = 0;
    while (kk < (posiFirst + len)) {
        array1[mm] = lastBeat[kk];
        mm = mm + 1;
        kk = kk + 1;
    }

    //% temp=abs(tempBeat(posSecond));
    k = floor(posSecond - errorLength);
    temp = fabsf(tempBeat[k]);
    i = k + 1;
    while (i < (posSecond + errorLength)) {
        if (fabsf(tempBeat[i]) > temp) {
            temp = fabsf(tempBeat[i]);
            k = i;
        }
        i = i + 1;
    }
    if (k != posSecond)
        posSecond = k;
    //end

    if (posSecond > len)
        kk = posSecond - len;
    else
        kk = 1;
    //end
    mm = 0;
    while (kk < (posSecond + len)) {
        array2[mm] = tempBeat[kk];
        mm = mm + 1;
        kk = kk + 1;
    }
    // [ correlation ] = correlationCoefficient( array1,array2 );
    float correlation = correlationCoefficient(array1, array2, len);

    return correlation;

}

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void templateCopy(const float *sourceBeat, float *destData) {

    int kk = FIDMARK - BEAT_MS60;
    int mm = 0;
    while (kk < FIDMARK + BEAT_MS90) {
        destData[mm] = sourceBeat[kk];
        mm = mm + 1;
        kk = kk + 1;
    }

}


//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void findKnots(float *y, int *knotsPosition, int *originalPosition, float *firstDiffer) {

    //utility tempPointer;
    float Xpc[sampleRate] = {0};
    // tempPointer.LowPassFilter(y,Xpc);
    LowPassFilter(y, Xpc);
    int lenth = sampleRate;
    float diff[sampleRate] = {0};//zeros(1,lenth);
    float absDiff[sampleRate] = {0};//zeros(1,lenth);

    int index = 0;
    while (index < (lenth - 1)) {
        diff[index] = Xpc[index + 1] - Xpc[index];
        absDiff[index] = fabsf(diff[index]);
        index = index + 1;
    }

    index = 0;
    float points[sampleRate] = {0};//=zeros(1,lenth);
    while (index < lenth) {
        points[index] = fabsf(diff[index] * Xpc[index]);//%array(index+1)-array(index);
        index = index + 1;
    }

    // %%%%%%%%%%%% find the smallest point  %%%%%%%%%%%%%%%%%%%
    int positionFirst = 0;
    int positionSecond = 0;
    index = BEAT_MS400;
    int currentMax = maxInteger;
    int currentMin = maxInteger;
    while (index > (BEAT_MS400 - BEAT_MS100)) {
        //%        temp=abs(array(index)/array(FIDMARK));
        if ((index < BEAT_MS360) && (absDiff[index] < (float) currentMax)) {
            currentMax = (int) absDiff[index];
            positionFirst = index;
        }

        if ((index < BEAT_MS360) && (points[index] < (float) currentMin)) {
            currentMin = (int) points[index];
            positionSecond = index;
        }
        index = index - 1;
    }


    *knotsPosition = 0;
    if ((positionFirst != 0) && (positionSecond != 0)) {
        if (fabsf(Xpc[positionFirst]) < (fabsf(Xpc[positionSecond]))) {
            *knotsPosition = positionFirst;
        } else {
            *knotsPosition = positionSecond;
        }
    }

    if (*knotsPosition > 1) {
        *firstDiffer = diff[*knotsPosition];
        // % originalPosition=y(knotsPosition);
        *originalPosition = (int) Xpc[*knotsPosition];
    }

}// end main

// int ecgMain::sinusIrregularityConut=0;
// int ecgMain::sinusArrestCount=0;
// int ecgMain::sinusBradycardiaCount=0;
// int ecgMain::sinusTachycardiaCount=0;
// int ecgMain::sickSinusSyndromeCount=0;
// int ecgMain::modeCycle=0;
// int ecgMain::sinoAtrialArray[10]={0};
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void sinusArrhythmia(sinusArrhythmia_CTX *ctx, int ppInterval, int *sinusArrhythmiaType, int *sinusIrregularityType,
                     int *sinusArrestType,
                     int *sinusBradycardiaType, int *sinusTachycardiaType, int *sickSinusSyndromeType, int init) {


    int *ctx_sinusIrregularityConut = &(ctx->sinusIrregularityConut);
    int *ctx_sinusArrestCount = &(ctx->sinusArrestCount);
    int *ctx_sinusBradycardiaCount = &(ctx->sinusBradycardiaCount);
    int *ctx_sinusTachycardiaCount = &(ctx->sinusTachycardiaCount);
    int *ctx_sickSinusSyndromeCount = &(ctx->sickSinusSyndromeCount);
    int *ctx_modeCycle = &(ctx->modeCycle);
    int *sinoAtrialArray = ctx->sinoAtrialArray;

    if (init) {
        (*ctx_sinusIrregularityConut) = 0;
        (*ctx_sinusArrestCount) = 0;
        (*ctx_sinusBradycardiaCount) = 0;
        (*ctx_sinusTachycardiaCount) = 0;
        (*ctx_sickSinusSyndromeCount) = 0;
        (*ctx_modeCycle) = 0;
        memset(sinoAtrialArray, 0, sizeof(int) * 10);
        return;
    }


    int arrayLength = 10;
    double delta = 0.2;
    (*ctx_modeCycle) = (*ctx_modeCycle) + 1;

    if (ppInterval < 0)
        return;
    //end
    int index;
    //if((*ctx_modeCycle)<=arrayLength)
    //{
    index = 0;
    while (index < (sinusArrLen - 1)) {
        sinoAtrialArray[index] = sinoAtrialArray[index + 1];
        index = index + 1;
    }
    sinoAtrialArray[sinusArrLen - 1] = ppInterval;
    //}
    // else
    if ((*ctx_modeCycle) >= sinusArrLen) {
        (*ctx_modeCycle) = 0;

        int sinoAtrialMaxPosition = 1;
        int sinoAtrialMax = sinoAtrialArray[0];
        int sinoAtrialMinPosition = 1;
        int sinoAtrialMin = sinoAtrialArray[0];

        index = 0;
        int diffCount = 0;
        int ppIntervalTotal = 0;
        int ppIntervalAverage = 0;
        while ((index < (sinusArrLen - 1))) {
            if (sinoAtrialArray[index] > sinoAtrialMax) {
                sinoAtrialMaxPosition = index;
                sinoAtrialMax = sinoAtrialArray[index];
            } else if (sinoAtrialArray[index] < sinoAtrialMin) {
                sinoAtrialMinPosition = index;
                sinoAtrialMin = sinoAtrialArray[index];
            }
            //end
            if (abs(sinoAtrialArray[index] - sinoAtrialArray[index + 1]) > BEAT_MS120) {
                diffCount = diffCount + 1;
            }
            ppIntervalTotal = ppIntervalTotal + sinoAtrialArray[index];
            index = index + 1;
        }


        ppIntervalAverage = ppIntervalTotal / sinusArrLen;
        index = sinoAtrialMaxPosition;
        while ((index < (sinusArrLen - 1)) && (index < sinoAtrialMinPosition)) {
            if (sinoAtrialArray[index] < (sinoAtrialArray[index + 1])) {
                int ppp = 0;
            } else
                break;
            //end
            index = index + 1;
        }

        if ((index < sinoAtrialMinPosition) && (sinoAtrialMax < 2 * sinoAtrialMin) && (sinoAtrialMax > BEAT_MS1200))
            *sinusArrhythmiaType = 1;
        //end

        // %%%%%%%%%%  modulus operation  %%%%%%%%%%%%%%%%%%
        int initial = 0;
        int tempSinoAtrialMax = sinoAtrialMax;
        int tempSinoAtrialMin = sinoAtrialMin;
        while ((tempSinoAtrialMax > tempSinoAtrialMin) && (tempSinoAtrialMin > 0)) {
            tempSinoAtrialMax = tempSinoAtrialMax - tempSinoAtrialMin;
            initial = initial + 1;
        }
        float surplus = (float) tempSinoAtrialMax / (float) sinoAtrialMin;
        if ((sinoAtrialMin > 0) && (sinoAtrialMax > BEAT_MS1200) && (surplus < delta))
            if ((initial >= 2) && (initial < 4))
                *sinusArrhythmiaType = 2;
        //end

        // %%%%%%%  sinus irregularity arrhythmia %%%%%%%%%%%%%%%%%
        if (sinoAtrialMin > 0 && sinusArrLen > 1 && (sinoAtrialMax - sinoAtrialMin) > BEAT_MS160)
            (*ctx_sinusIrregularityConut) = (*ctx_sinusIrregularityConut) + 1;

        // %%%%%%%%%%%% sinus arrest  %%%%%%%%%%%%%%%%%
        if ((sinoAtrialMin > 0) && (sinoAtrialMax > BEAT_MS2000))
            (*ctx_sinusArrestCount) = (*ctx_sinusArrestCount) + 1;

        //%%%%%%%%%%%%%%%%%%%    Sinus Bradycardia      %%%%%%%%%%%%%%%%%%%
        if (ppIntervalAverage > BEAT_MS1100)
            (*ctx_sinusBradycardiaCount) = (*ctx_sinusBradycardiaCount) + 1;
        // end
        // %%%%%%%%%%%%%%%%%%% Sinus Tachycardia  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sick Sinus Syndrome
        if (ppIntervalAverage < BEAT_MS600)
            (*ctx_sinusTachycardiaCount) = (*ctx_sinusTachycardiaCount) + 1;

        // %%%%%%%%%%%%%%%%%%% Sick Sinus Syndrome  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((sinusBradycardiaType) && (sinusArrestType))
            (*ctx_sickSinusSyndromeCount) = (*ctx_sickSinusSyndromeCount) + 1;

    }

    *sinusArrestType = (*ctx_sinusArrestCount);
    *sinusIrregularityType = (*ctx_sinusIrregularityConut);
    *sinusBradycardiaType = (*ctx_sinusBradycardiaCount);
    *sinusTachycardiaType = (*ctx_sinusTachycardiaCount);
    *sickSinusSyndromeType = (*ctx_sickSinusSyndromeCount);

}

//
// int ecgMain:: vpcString[vTachyLength]={0};
// int ecgMain:: vpcBigeminy=0;
// int ecgMain:: vpcTrigeminy=0;
// int ecgMain:: vpcBigString[vTachyLength+1]={0};
// int ecgMain:: vpcTriString[vTachyLength*2]={0};
// int ecgMain:: longestPvcTachycardiaCount=0;
// int ecgMain:: pvcTachycardiaSegment=0;
// int ecgMain:: pvcPointerCount=0;
// int ecgMain:: pvcBiSegmentCount=0;
// int ecgMain:: pvcBiLongest=0;
// int ecgMain:: fastestPvcTachycardia=0;
// int ecgMain:: doublePvcCount=0;
// int ecgMain:: doubleArray[vTachyLength]={0};
// int ecgMain:: rrSum=0;
// int ecgMain:: rrString[vTachyLength]={0};
// int ecgMain:: biFlag=0;
// int ecgMain:: flag=0;
// int ecgMain:: count=0;
// int ecgMain:: triCount=0;
// int ecgMain:: triFisrFlag=0;
// int ecgMain:: triSecongFlag=0;
// int ecgMain:: pvcTriSegmentCount=0;
// int ecgMain:: pvcTriLongest=0;
// int ecgMain:: AVArray[10]={0};
// int ecgMain:: firstDegreeBlockCount=0;

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int firstDegreeAVBlock(firstDegreeAVBlock_CTX *ctx, int prInterval, int init) {

    int *AVArray = ctx->AVArray;
    int *ctx_firstDegreeBlockCount = &(ctx->firstDegreeBlockCount);


    if (init) {
        memset(AVArray, 0, sizeof(int) * 10);
        (*ctx_firstDegreeBlockCount) = 0;

        return 0;
    }

    int arrayLength = 10;
    int index = 0;
    int prIntervalSum = 0;
    while ((index < arrayLength - 1)) {
        AVArray[index] = AVArray[index + 1];
        index = index + 1;
    }
    AVArray[arrayLength - 1] = prInterval;

    index = 0;
    int number = 0;
    while ((index < arrayLength)) {
        prIntervalSum = prIntervalSum + AVArray[index];
        if ((((index + 2) <= arrayLength) && AVArray[index] > BEAT_MS200) && (AVArray[index + 2] > BEAT_MS200)) {
            number = number + 1;
        }
        //% if((AVArray(index)<1)||(AVArray(index)<BEAT_MS200))
        if (AVArray[index] < 1) {
            prIntervalSum = 0;
            break;
        }
        index = index + 1;

    }

    if (((index >= arrayLength) && ((float) prIntervalSum / (float) arrayLength) > BEAT_MS230) || (number >= 5))
        (*ctx_firstDegreeBlockCount) = (*ctx_firstDegreeBlockCount) + 1;

//     if((*ctx_firstDegreeBlockCount)>3000)
//     {
//       (*ctx_firstDegreeBlockCount)=3000;

//     }
    //%     AVBlockFlag=(*ctx_firstDegreeBlockCount);
    index = 0;
    while (index < arrayLength) {
        AVArray[index] = 0;//zeros(1,arrayLength);
        index++;
    }
    //AVBlockFlag=(*ctx_firstDegreeBlockCount);
    return (*ctx_firstDegreeBlockCount);
    // end

}


//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void
ST_Morphology(const float *array, int begin, int finish, float amplifierGainConst, int len, int tFirstPeak, float Crit,
              int *type) {
    // #define  len  80

    //%%%%%%%%%%%%% ST Segment Analysis %%%%%%%%%%%%%%%%%%%%
    *type = 0;
    float heightDiffer[80] = {0};//zeros(1,len);
    float y[80] = {0};//zeros(1,len);
    // p=zeros(1,len);

    float max = array[begin + 1];
    float min = array[begin + 1];
    int maxPosition = begin;
    int minPosition = begin;
    int positiveNum = 0;
    int negativeNum = 0;
    int k = 1;
    int i = begin;

    while (i <= finish) {
        // %%%%  y  straight line   %%%%%
        y[i] = array[begin] +
               (float) (i - begin) * ((float) (array[tFirstPeak] - array[begin])) / (float) (tFirstPeak - begin);
        heightDiffer[k] = array[i] - y[i];
        if (heightDiffer[k] > 0)
            positiveNum = positiveNum + 1;
        else
            negativeNum = negativeNum + 1;
        // end
        if (heightDiffer[k] > max) {
            max = heightDiffer[k];
            maxPosition = k;
        } else {
            if (heightDiffer[k] < min)
                min = heightDiffer[k];
            minPosition = k;
        }// end
        i = i + 1;
        k = k + 1;
    }// end

    //%%%%%%%%  straight line  ST Type %%%%%%%%%%%%%%%%%%%%
    float stSlope = 0;
    if (begin != finish) {
        stSlope = (float) 500 * (array[finish] - array[begin]) / (float) (finish - begin);
    }
    if (stSlope > 10 * Crit)
        *type = 1;      // %%  type=1;  upslope
    else if ((stSlope) > 10 * Crit)
        *type = 2;   // %%  type=2;  downslope
    else
        *type = 3; // %%  type=3;  straight line
    // end
    //en
    // %%%%%%%%  curve  ST Type %%%%%%%%%%%%%%%%%%%%
    float positivePossibility = (float) positiveNum / (float) (len + 1);
    float negativePossibility = (float) negativeNum / (float) (len + 1);

    if (fabsf(max - min) < Crit / 3)
        *type = 4;        // %%strainghtLine=1;
    else if (positivePossibility > 0.7)
        *type = 5;     // %% ST_Convext=1;
    else if (negativePossibility > 0.7)
        *type = 6;   // %%  ST_Concave=1;
    //end
    // end
    // end

}


//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ST_SegmentAnalysis(float *array, int rr, int pStart, int pOffset, int onset, int offset, float amplifierGainConst,
                        int tFirstPeak, float Crit, float leadStandard, int *ST_MorphologyType, int *ST_Change,
                        float *ST_ElavatorAmph, int *type) {
    *ST_Change = 0;

    //[knotsPosition,originalPosition,firstDiffer] = findKnots(array);
    //int len=500;
    float y[Fs] = {0};
    arrayCopy(array, y, Fs);
    // % %%%% delete the p wave %%%%%%

    int i = pStart;
    while (((i <= pOffset) && (pStart != 0))) {
        y[i] = y[pStart];
        i = i + 1;
    }//end

    // %%%%%%%% find the electric level %%%%%%%%%%%%%
    i = 0;
    float pointTotal = 0;
    float isoelectricValue;
    while (i <= onset) {
        pointTotal = pointTotal + y[i];
        i = i + 1;
    }//end
    float iosTemp = (float) pointTotal / (float) onset;

    if (fabsf(y[onset]) < fabsf(iosTemp)) {
        isoelectricValue = y[onset];
    } else {
        isoelectricValue = iosTemp;
    }//end
    // if(abs(originalPosition)<abs(y(onset)))

    //   isoelectricValue=abs(originalPosition);
    //  end


    float RRav = (float) rr / Fs;
    float HR = (float) 60 / RRav;
    int ST_Offset = 0;
    if (HR >= 140)
        ST_Offset = BEAT_MS40;
    else if (HR >= 120)
        ST_Offset = BEAT_MS60;
    else if (HR >= 110)
        ST_Offset = BEAT_MS64;
    else if (HR >= 100)
        ST_Offset = BEAT_MS72;
    else if (HR < 100)
        ST_Offset = BEAT_MS80;
    //end

    int ST_MeasuringPoint = offset + ST_Offset;
    int len = ST_MeasuringPoint - offset;
    //  %array(onset)=0;
    *ST_ElavatorAmph = (float) (y[ST_MeasuringPoint] - isoelectricValue) / amplifierGainConst;
    if (*ST_ElavatorAmph > leadStandard) {
        *ST_Change = (int) *ST_ElavatorAmph;
    }//end
    int mType = 0;
    if ((len > 1) && (offset != ST_MeasuringPoint)) {

        ST_Morphology(array, offset, ST_MeasuringPoint, amplifierGainConst, len, tFirstPeak, Crit, &mType);
        *ST_MorphologyType = mType;

    }//end

    *type = mType;

}


//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ventricularDataCopy(int doublePvcCount, int pvcBiSegmentCount, int pvcTriSegmentCount, int pvcTachycardiaSegment,
                         int longestPvcTachycardiaCount, int fastestPvcTachycardia, int *doublePvc,
                         int *vpcBigeminyCount, int *vpcTrigeminyCount, int *ventricularTachycardiaConut,
                         int *longestPvcTachycardia, int *fastestPvcTachy) {
    *doublePvc = doublePvcCount;
    *vpcBigeminyCount = pvcBiSegmentCount;
    *vpcTrigeminyCount = pvcTriSegmentCount;
    *ventricularTachycardiaConut = pvcTachycardiaSegment;
    *longestPvcTachycardia = longestPvcTachycardiaCount;
    *fastestPvcTachy = fastestPvcTachycardia;

}

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void clearBuffer(int *vpcTriString, int *doubleArray, int *vpcBigString) {
    int ii = 0;
    while (ii < 6) {
        vpcTriString[ii] = 0;
        ii = ii + 1;
    }//end

    ii = 0;
    while (ii < 3) {
        doubleArray[ii] = 0;
        ii = ii + 1;
    }

    ii = 0;
    while (ii < 4) {
        vpcBigString[ii] = 0;
        ii = ii + 1;
    }

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ventricularPrematureCount(ventricularPrematureCount_CTX *ctx, int vpcFlag, int rr, int *doublePvc,
                               int *vpcBigeminyCount, int *vpcTrigeminyCount,
                               int *ventricularTachycardiaCount, int *longestPvcTachycardia, int *fastestPvcTachy,
                               int init) {
    //  %%%%%%% ventricular tachycardia (*ctx_count) %%%%%%

    int *vpcString = ctx->vpcString; // static int vpcString[vTachyLength] = {0};
    int *vpcBigString = ctx->vpcBigString; // static int vpcBigString[vTachyLength + 1] = {0};
    int *vpcTriString = ctx->vpcTriString; // static int vpcTriString[vTachyLength * 2] = {0};
    int *doubleArray = ctx->doubleArray; // static int doubleArray[vTachyLength] = {0};
    int *rrString = ctx->rrString; // static int rrString[vTachyLength] = {0};

    int *ctx_vpcBigeminy = &(ctx->vpcBigeminy);
    int *ctx_vpcTrigeminy = &(ctx->vpcTrigeminy);
    int *ctx_longestPvcTachycardiaCount = &(ctx->longestPvcTachycardiaCount);
    int *ctx_pvcTachycardiaSegment = &(ctx->pvcTachycardiaSegment);
    int *ctx_pvcPointerCount = &(ctx->pvcPointerCount);
    int *ctx_pvcBiSegmentCount = &(ctx->pvcBiSegmentCount);
    int *ctx_pvcBiLongest = &(ctx->pvcBiLongest);
    int *ctx_fastestPvcTachycardia = &(ctx->fastestPvcTachycardia);
    int *ctx_doublePvcCount = &(ctx->doublePvcCount);
    int *ctx_rrSum = &(ctx->rrSum);
    int *ctx_biFlag = &(ctx->biFlag);
    int *ctx_flag = &(ctx->flag);
    int *ctx_count = &(ctx->count);
    int *ctx_triCount = &(ctx->triCount);
    int *ctx_triFisrFlag = &(ctx->triFisrFlag);
    int *ctx_triSecongFlag = &(ctx->triSecongFlag);
    int *ctx_pvcTriSegmentCount = &(ctx->pvcTriSegmentCount);
    int *ctx_pvcTriLongest = &(ctx->pvcTriLongest);

    if (init) {

        // vpcString[vTachyLength]={0};
        memset(vpcString, 0, sizeof(int) * vTachyLength);
        (*ctx_vpcBigeminy) = 0;
        (*ctx_vpcTrigeminy) = 0;
        //vpcBigString[vTachyLength+1]={0};
        memset(vpcBigString, 0, sizeof(int) * (vTachyLength + 1));
        // vpcTriString[vTachyLength*2]={0};
        memset(vpcTriString, 0, sizeof(int) * (vTachyLength * 2));
        (*ctx_longestPvcTachycardiaCount) = 0;
        (*ctx_pvcTachycardiaSegment) = 0;
        (*ctx_pvcPointerCount) = 0;
        (*ctx_pvcBiSegmentCount) = 0;
        (*ctx_pvcBiLongest) = 0;
        (*ctx_fastestPvcTachycardia) = 0;
        (*ctx_doublePvcCount) = 0;
        // doubleArray[vTachyLength]={0};
        memset(doubleArray, 0, sizeof(int) * (vTachyLength));
        (*ctx_rrSum) = 0;
        // rrString[vTachyLength]={0};
        memset(rrString, 0, sizeof(int) * (vTachyLength));
        (*ctx_biFlag) = 0;
        (*ctx_flag) = 0;
        (*ctx_count) = 0;
        (*ctx_triCount) = 0;
        (*ctx_triFisrFlag) = 0;
        (*ctx_triSecongFlag) = 0;
        (*ctx_pvcTriSegmentCount) = 0;
        (*ctx_pvcTriLongest) = 0;

        return;

    }


    float temp;
    int ii;
    if ((*ctx_pvcPointerCount) == 0) {
        ii = 0;
        while (ii < (vTachyLength - 1)) {
            vpcString[ii] = vpcString[ii + 1];
            rrString[ii] = rrString[ii + 1];
            ii = ii + 1;
        }
        if (vpcFlag == 1) {
            vpcString[vTachyLength - 1] = 1;
            rrString[vTachyLength - 1] = rr;
        } else {
            vpcString[vTachyLength - 1] = 0;
            rrString[vTachyLength - 1] = 0;
        }
    } else {
        if ((*ctx_pvcPointerCount) > (*ctx_longestPvcTachycardiaCount)) {
            (*ctx_longestPvcTachycardiaCount) = (*ctx_pvcPointerCount);
        }
        if (vpcFlag == 1) {
            // %%%%%% long tachycardia %%%%%%
            (*ctx_flag) = 1;
            (*ctx_pvcPointerCount) = (*ctx_pvcPointerCount) + 1;
            (*ctx_rrSum) = (*ctx_rrSum) + rr;
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;
        } else {
            temp = (float) 2 * (float) (*ctx_rrSum) / (float) (*ctx_pvcPointerCount);
            float period = (float) 60000 / temp;
            if (period > (float) (*ctx_fastestPvcTachycardia)) {
                (*ctx_fastestPvcTachycardia) = (int) period;
            }
            (*ctx_pvcPointerCount) = 0;
            (*ctx_pvcTachycardiaSegment) = (*ctx_pvcTachycardiaSegment) + 1;
            (*ctx_flag) = 0;
            ii = 0;
            while (ii < 3) {
                vpcString[ii] = 0;
                rrString[ii] = 0;
                ii = ii + 1;
            }
            clearBuffer(vpcTriString, doubleArray, vpcBigString);
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;
        }
    }

    if (((*ctx_flag) == 0) && (vpcString[0] == 1) && (vpcString[1] == 1) && (vpcString[2] == 1)) {
        (*ctx_pvcPointerCount) = 3;
        (*ctx_rrSum) = rrString[0] + rrString[1] + rrString[2];
        // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
        ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                            (*ctx_pvcTachycardiaSegment),
                            (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                            vpcBigeminyCount,
                            vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
        return;
    }//end

    // %%%%%%% ventricular tachycardia (*ctx_count)  end %%%%%%;

    //   %%%%%%%  bigeminy counter  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

    if ((*ctx_vpcBigeminy) == 0) {
        ii = 0;
        while (ii < 3) {
            vpcBigString[ii] = vpcBigString[ii + 1];
            ii = ii + 1;
        }//end
        if (vpcFlag == 1) {
            vpcBigString[3] = 1;
        } else {
            vpcBigString[3] = 0;
        }//end
    } else {
        (*ctx_count) = (*ctx_count) + 1;
        if (((*ctx_count) == 1) && (vpcFlag == 0)) {
            (*ctx_biFlag) = 1;
            //[doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;

        }

        if ((vpcFlag == 1) && ((*ctx_biFlag) == 1) && ((*ctx_count) == 2)) {
            (*ctx_count) = 0;
            (*ctx_biFlag) = 0;
            (*ctx_vpcBigeminy) = (*ctx_vpcBigeminy) + 1;
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;
        } else {
            (*ctx_pvcBiSegmentCount) = (*ctx_pvcBiSegmentCount) + 1;
            if ((*ctx_pvcBiLongest) < (*ctx_vpcBigeminy)) {
                (*ctx_pvcBiLongest) = (*ctx_vpcBigeminy);
            }
            (*ctx_vpcBigeminy) = 0;
            (*ctx_count) = 0;

            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;
        }
    }

    if ((vpcBigString[0] == 0) && (vpcBigString[1] == 1) && (vpcBigString[2] == 0) && (vpcBigString[3] == 1)) {
        (*ctx_vpcBigeminy) = (*ctx_vpcBigeminy) + 1;
        ii = 0;
        while (ii < 4) {
            vpcBigString[ii] = 0;
            ii = ii + 1;
        }//end
        (*ctx_biFlag) = 0;
        (*ctx_count) = 0;
        // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
        ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                            (*ctx_pvcTachycardiaSegment),
                            (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                            vpcBigeminyCount,
                            vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
        return;
    }//end

    //%%%%%%% bigeminy counter  end  %%%

    ii = 0;
    while (ii < 2) {
        doubleArray[ii] = doubleArray[ii + 1];
        ii = ii + 1;
    }
    if (vpcFlag == 1)
        doubleArray[2] = 1;
    else {
        doubleArray[2] = 0;
    }//end

    //    %%%%%%%%%%% double pvc %%%%%%%%%%%%
    if ((doubleArray[0] == 1) && (doubleArray[1] == 1)) {
        (*ctx_doublePvcCount) = (*ctx_doublePvcCount) + 1;
        ii = 0;
        while (ii < 3) {
            doubleArray[ii] = 0;
            ii = ii + 1;
        }//end
        // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
        ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                            (*ctx_pvcTachycardiaSegment),
                            (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                            vpcBigeminyCount,
                            vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
        return;
    }//end


    //   %%%%%%%%%%%%%%% Trigeminy %%%%%%%%%%%
    if ((*ctx_vpcTrigeminy) == 0) {
        ii = 0;
        while (ii < 5) {
            vpcTriString[ii] = vpcTriString[ii + 1];
            ii = ii + 1;
        }//end
        if (vpcFlag == 1)
            vpcTriString[5] = 1;
        else {
            vpcTriString[5] = 0;
        }//end
    } else {
        (*ctx_triCount) = (*ctx_triCount) + 1;
        if (((*ctx_triCount) == 1) && (vpcFlag == 0)) {
            (*ctx_triFisrFlag) = 1;
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;
        } else {
            if (((*ctx_triCount) == 2) && (vpcFlag == 0)) {
                (*ctx_triSecongFlag) = 1;
                //[doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
                ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                    (*ctx_pvcTachycardiaSegment),
                                    (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                    vpcBigeminyCount,
                                    vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia,
                                    fastestPvcTachy);
                return;
            }//end
        }// end

        if (((*ctx_triFisrFlag) == 1) && ((*ctx_triSecongFlag) == 1) && (vpcFlag == 1) && ((*ctx_triCount) == 3)) {
            (*ctx_vpcTrigeminy) = (*ctx_vpcTrigeminy) + 1;
            (*ctx_triFisrFlag) = 0;
            (*ctx_triSecongFlag) = 0;
            (*ctx_triCount) = 0;
        } else {
            (*ctx_pvcTriSegmentCount) = (*ctx_pvcTriSegmentCount) + 1;
            if ((*ctx_pvcTriLongest) < (*ctx_vpcTrigeminy)) {
                (*ctx_pvcTriLongest) = (*ctx_vpcTrigeminy);
            }//end
            (*ctx_vpcTrigeminy) = 0;
            (*ctx_triCount) = 0;
            //[doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
            ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                                (*ctx_pvcTachycardiaSegment),
                                (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc,
                                vpcBigeminyCount,
                                vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
            return;

        }//end


    }// end

    if ((vpcTriString[0] == 0) && (vpcTriString[1] == 0) && (vpcTriString[2] == 1) && (vpcTriString[3] == 0) &&
        (vpcTriString[4] == 0) && (vpcTriString[5] == 1)) {
        (*ctx_vpcTrigeminy) = (*ctx_vpcTrigeminy) + 1;
        ii = 0;
        while (ii < 6) {
            vpcTriString[ii] = 0;
            ii = ii + 1;
        }//end
        (*ctx_triCount) = 0;
    }//end
    //[doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( (*ctx_doublePvcCount)  ,(*ctx_pvcBiSegmentCount),(*ctx_pvcTriSegmentCount),(*ctx_pvcTachycardiaSegment),(*ctx_longestPvcTachycardiaCount),(*ctx_fastestPvcTachycardia));
    ventricularDataCopy((*ctx_doublePvcCount), (*ctx_pvcBiSegmentCount), (*ctx_pvcTriSegmentCount),
                        (*ctx_pvcTachycardiaSegment),
                        (*ctx_longestPvcTachycardiaCount), (*ctx_fastestPvcTachycardia), doublePvc, vpcBigeminyCount,
                        vpcTrigeminyCount, ventricularTachycardiaCount, longestPvcTachycardia, fastestPvcTachy);
    //  %%%%%%%%%%%%%%% Trigeminy  end  %%%%%%%%%%%

}

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float templateSeek(PvcLibInfor *array, float *source) {

    float temp[100] = {0};//=zeros(1,100);
    float maxCoef = 0;
    int len = array->queueLength;

    int i = 0, k;
    while (i < len) {
        k = 0;
        while (k < 100) {
            temp[k] = array->informArray[i][k];
            k = k + 1;
        }
        float correlation = correlationCoefficient(temp, source, 100);
        if (correlation > maxCoef) {
            maxCoef = correlation;
        }

        i = i + 1;

    }//end

    return maxCoef;

}

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float shannonDistribution(float *rrArray, int len) {
#define  shanLen  30

    //  BEAT_MS1200=600;
    //  BEAT_MS120=60;
    int deltaRRmax = BEAT_MS1200;
    int deltaRRmin = -BEAT_MS1200;
    int RRImaxThreshold = BEAT_MS120;
    int RRImaxCount = 0;
    int histNumber = 16;
    //    histWidth=16;
    int tailHeadLen = 2;//%5;%8;

    float tempArray[shanLen] = {0};
    float templateArray[shanLen] = {0};

    arrayCopy(rrArray, tempArray, len);
    arrayCopy(rrArray, templateArray, len);
    DataSort(templateArray, len);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int jj = 0;
    while ((jj < len) && (templateArray[jj] == 0)) {
        jj = jj + 1;
    }

    int kk = jj;
    int ll = 0;
    int ii;
    float tempData;
    while ((ll < tailHeadLen) && (kk < len)) {
        tempData = templateArray[kk];
        ii = 0;
        while (ii < len) {
            if (tempData == tempArray[ii] && (tempArray[ii] != 0)) {
                tempArray[ii] = 0;
                ii = len;
            }
            ii = ii + 1;
        }
        kk = kk + 1;
        ll = ll + 1;
    }
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kk = len;
    ll = 0;
    while (ll < tailHeadLen) {
        tempData = templateArray[kk];
        ii = 0;
        while (ii < len) {
            if (tempData == tempArray[ii] && (tempArray[ii] != 0)) {
                tempArray[ii] = 0;
                ii = len;
            }
            ii = ii + 1;
        }
        kk = kk - 1;
        ll = ll + 1;
    }//end
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int effectiveLength = 0;
    float rrIntervalArray[shanLen] = {0};//zeros(1,len);
    ii = 0;
    while (ii < len) {
        if (tempArray[ii] != 0) {
            rrIntervalArray[effectiveLength] = tempArray[ii];
            effectiveLength = effectiveLength + 1;
        }
        ii = ii + 1;
    }
    // rrArray=rrIntervalArray;
    arrayCopy(rrIntervalArray, templateArray, shanLen);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int nbins, length;
    float diffRRarray[shanLen] = {0};
    if (effectiveLength < 4) {
        float shannonEntropy = 0;
        return shannonEntropy;
    } else {
        nbins = effectiveLength;//%histWidth ;% length;%% 40;
        length = effectiveLength;
        // float diffRRarray[shannLen]={0};//zeros(1,nbins-1);
    }//end

    int index = 0;
    while (index < (length - 1)) {
        diffRRarray[index] = templateArray[index + 1] - templateArray[index];
        index = index + 1;
    }
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    index = 0;
    float diffSum = 0;
    while (index < nbins) {
        diffSum = diffSum + diffRRarray[index];
        index = index + 1;
    }//end

    index = 0;
    float diffMean = (float) diffSum / (float) nbins;
    float diffSquareSum = 0;
    float RMSSD = 0;
    while (index < nbins) {
        diffSquareSum = diffSquareSum + (diffRRarray[index] - diffMean) * (diffRRarray[index] - diffMean);
        if (fabsf(diffRRarray[index] - diffMean) > (float) (RRImaxThreshold)) {
            RRImaxCount = RRImaxCount + 1;
        }
        index = index + 1;
    }
    float sumRR = 0;
    index = 0;
    while (index < (nbins - 1)) {
        sumRR = sumRR + templateArray[index];
//     int temp1=rrArray[index+1];
//     int temp2=rrArray[index];
//     int temp3=(rrArray[index+1]-rrArray[index]);

        RMSSD = RMSSD +
                (templateArray[index + 1] - templateArray[index]) * (templateArray[index + 1] - templateArray[index]);
        index = index + 1;
    }
    float rrMean = sumRR / (float) nbins;
    RMSSD = sqrtf(RMSSD / (float) (nbins - 1));
    float tempRR = RMSSD;
    RMSSD = RMSSD / rrMean;
//  %%%%%%%%%%%%%  build hist %%%%%%%%%%%%%%%%%%%%%%
    int histPropability[histWidth] = {0};//zeros(1,histWidth);
    index = 0;
    int num;
    float temp;
    int diffNum = length;
    while (index < histNumber) {
        num = 0;
        while (num < diffNum) {
            temp = floorf(
                    (float) histNumber * (diffRRarray[num] - (float) deltaRRmin) / (float) (deltaRRmax - deltaRRmin));
            if (temp == (float) index) {
                histPropability[index] = histPropability[index] + 1;
            }
            num = num + 1;
        }
        index = index + 1;
    }//end
    float Shannon1 = 0;
    index = 0;
    float p;
    int histValidatedCount = 0;
    int histMaxCount = histPropability[index];
    while (index < histWidth) {
        if (histPropability[index] == 0) {
            p = 1;
        } else {
            p = (float) histPropability[index] / (float) effectiveLength;
            histValidatedCount = histValidatedCount + 1;
        }
        if (histPropability[index] > histMaxCount) {
            histMaxCount = histPropability[index];
        }
        Shannon1 = Shannon1 + (float) (p * log2f(p)) / log2f(1 / (float) histNumber);
        index = index + 1;
    }//end
    float shannonEntropy = RMSSD;

    return shannonEntropy;

}


//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float correlationCoefficient(const float *array1, const float *array2, int len) {
    float sumX = 0;
    float sumY = 0;
    float sumXX = 0;
    float sumXY = 0;
    float sumYY = 0;

    int k = 0;
    while (k < len) {
        sumX = sumX + array1[k];
        sumY = sumY + array2[k];
        sumXX = sumXX + array1[k] * array1[k];
        sumXY = sumXY + array1[k] * array2[k];
        sumYY = sumYY + array2[k] * array2[k];
        k = k + 1;
    }

    float correlation =
            ((float) len * sumXY - sumX * sumY) /
            (sqrtf((float) len * sumXX - sumX * sumX) * sqrtf((float) len * sumYY - sumY * sumY));

    return correlation;

}


//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float FuzzyReasoning(int classification, float inputData) {
    float returnValue = 0;
    float constA = 0, constB = 0, constC = 0, constD = 0;

    if (classification == 1) {
        constA = 20;//%-0.4;
        constB = 58;//%%75;%%
        constC = 120;//%90;%85;
        constD = 180;
    } else if (classification == 2) {
        constA = 0;
        constB = 140;
        constC = 210;//%165;
        constD = 250;
    } else if (classification == 3) {
        constA = -0.8f;
        constB = -0.65f;//%-0.65;%-0.4;
        constC = -0.196f;//%%-0.20;
        constD = 2;
    } else if (classification == 4) {
        constA = -0.6f;
        constB = 0.1f;
        constC = 0.5f;
        constD = 2.50f;
    } else if (classification == 5) {
        constA = 0.0f;
        constB = 0.7f;
        constC = 0.8f;
        constD = 2.50f;
    } else if (classification == 6) {
        constA = 0;
        constB = 1.75f;//%2;
        constC = 4;
        constD = 6;
    }


    float returnValueA = 0;
    if (inputData <= constA) {
        returnValueA = 0;
    } else if ((inputData > constA) && (inputData <= (float) (constA + constB) / 2)) {
        returnValueA = 2 * (inputData - constA) * (inputData - constA) / ((constB - constA) * (constB - constA));
    } else if ((inputData >= (float) (constB + constA) / 2) && (inputData < constB)) {
        returnValueA = 1 - 2 * (constB - inputData) * (constB - inputData) / ((constB - constA) * (constB - constA));
    } else if (inputData >= constB) {
        returnValueA = 1;
    }//end
    // end
    // end
    // end

    float returnValueB = 0;
    if (inputData <= constC) {
        returnValueB = 1;
    } else {
        if ((inputData > constC) && (inputData <= (float) (constC + constD) / 2)) {
            if (classification == 3) {
                if ((inputData > -0.2) && (inputData < -0.15)) {
                    returnValueB = 0.75f;
                } else {
                    returnValueB =
                            2 * (constC - inputData) * (constC - inputData) / ((constD - constC) * (constD - constC));
                }//end
            } else {
                returnValueB =
                        1 - 2 * (constC - inputData) * (constC - inputData) / ((constD - constC) * (constD - constC));
            }// end
        } else {
            if ((inputData >= (float) (constC + constD) / 2) && (inputData < constD)) {
                if (classification == 3) {
                    returnValueB =
                            0.5f * (inputData - constD) * (inputData - constD) /
                            ((constD - constC) * (constD - constC));
                } else {
                    returnValueB =
                            2 * (inputData - constD) * (inputData - constD) / ((constD - constC) * (constD - constC));
                }
            } else {
                if (inputData >= constD) {
                    returnValueB = 0;
                }//end
            }//end
        }//  end
    }// end

    returnValue = returnValueA * returnValueB;
    return returnValue;

}

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
float templateCompute(const float *array, int number) {
    float AverResult = 0;
    int i = 0;
    while (i < number) {
        AverResult = AverResult + array[i];
        i = i + 1;
    }
    AverResult = (float) AverResult / (float) number;
    return AverResult;

}


//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void DataSort(float *dataArray, int length) {
    int index_i, index_j;
    float temp;
    // for index_i=1:length-1
    for (index_i = 0; index_i < length - 1; index_i++) {
        // for index_j=index_i+1:length
        for (index_j = index_i + 1; index_j < length; index_j++) {
            if (dataArray[index_j] < dataArray[index_i]) {
                temp = dataArray[index_i];
                dataArray[index_i] = dataArray[index_j];
                dataArray[index_j] = temp;
            }//end
        }//end
    }//end

}

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void arrayCopy(const float *sourceArray, float *destArray, int len) {
    int ii = 0;
    while (ii < len) {
        destArray[ii] = sourceArray[ii];
        ii++;
    }
}

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void
shannonEntropyCompute(float *rrArray, const float *prDistribution, float *shannonEntropy, float *RMSSD,
                      float *prIntervalDiff,
                      float *RRImaxCount) {

    int deltaRRmax = BEAT_MS1200;
    int deltaRRmin = -BEAT_MS1200;
    int RRImaxThreshold = BEAT_MS120;
    *RRImaxCount = 0;
    int histNumber = 16;
    //   int histWidth=16;
    int len = 64;//%128;
    int tailHeadLen = 5;//%5


    float tempArray[shannLen] = {0};
    //tempArray=rrArray;
    arrayCopy(rrArray, tempArray, len);

    DataSort(rrArray, len);

    int jj = 0;
    while ((jj < len) && (rrArray[jj] == 0)) {
        jj = jj + 1;
    }

    float tempData;
    int kk = jj;
    int ll = 0;
    int ii;
    while ((ll < tailHeadLen) && (kk < len)) {
        tempData = rrArray[kk];
        ii = 0;
        while (ii < len) {
            if (tempData == tempArray[ii] && (tempArray[ii] != 0)) {
                tempArray[ii] = 0;
                ii = len;
            }
            ii = ii + 1;
        }
        kk = kk + 1;
        ll = ll + 1;
    }//end


    kk = len - 1;
    ll = 0;
    while (ll < tailHeadLen) {
        tempData = rrArray[kk];
        ii = 0;
        while (ii < len) {
            if (tempData == tempArray[ii] && (tempArray[ii] != 0)) {
                tempArray[ii] = 0;
                ii = len;
            }
            ii = ii + 1;
        }
        kk = kk - 1;
        ll = ll + 1;
    }//end

    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int effectiveLength = 0;
    float rrIntervalArray[shannLen] = {0};//zeros(1,len);
    ii = 0;
    while (ii < len) {
        if (tempArray[ii] != 0) {
            rrIntervalArray[effectiveLength] = tempArray[ii];
            effectiveLength = effectiveLength + 1;
        }
        ii = ii + 1;
    }//end

    // rrArray=rrIntervalArray;
    arrayCopy(rrIntervalArray, rrArray, shannLen);
    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (effectiveLength <= (tailHeadLen + 1)) {
        shannonEntropy = 0;
        *RMSSD = 0;
        prIntervalDiff = 0;
        RRImaxCount = 0;
        return;
    }

    int nbins = effectiveLength;//%histWidth ;//% length;%% 40;
    int length = effectiveLength;
    float diffRRarray[shannLen] = {0};//zeros(1,nbins-1);

    int index = 0;
    while (index < (length - 1)) {
        diffRRarray[index] = rrArray[index + 1] - rrArray[index];
        index = index + 1;
    }//end
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    index = 0;
    int prIntervaLength = 0;
    float prIntervalSum = 0;
    while (index < len) {
        prIntervalSum = prIntervalSum + prDistribution[index];
        if (prDistribution[index] != 0) {
            prIntervaLength = prIntervaLength + 1;
        }
        index = index + 1;
    }//end

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float prArray[shannLen] = {0};//zeros(1,prIntervaLength);
    index = 0;
    ii = 0;
    while (index < len) {
        if (prDistribution[index] != 0) {
            prArray[ii] = prDistribution[index];
            ii = ii + 1;
        }
        index = index + 1;
    }//end
    float prIntervalMean = (float) prIntervalSum / (float) prIntervaLength;


//  $$$$$$$$$$$$$$$$$$$

    float prIntervalStandardDiff = 0;
    index = 0;
    while (index < prIntervaLength) {
        prIntervalStandardDiff = prIntervalStandardDiff +
                                 (prDistribution[index] - prIntervalMean) * (prDistribution[index] - prIntervalMean);
        index = index + 1;
    }
    *prIntervalDiff = sqrtf(prIntervalStandardDiff / (float) prIntervaLength);

//  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    index = 0;
    float diffSum = 0;
    // %prIntervalSum=0;
    while (index < nbins) {
        diffSum = diffSum + diffRRarray[index];
        // %  prIntervalSum=prIntervalSum+prDistribution(index);
        index = index + 1;
    }//end

    index = 0;
    float diffMean = diffSum / (float) nbins;
    // %prIntervalMean=prIntervalSum/nbins;
    prIntervalStandardDiff = 0;
    float diffSquareSum = 0;
    *RMSSD = 0;
    while (index < nbins) {
        diffSquareSum = diffSquareSum + (diffRRarray[index] - diffMean) * (diffRRarray[index] - diffMean);
        prIntervalStandardDiff = prIntervalStandardDiff +
                                 (prDistribution[index] - prIntervalMean) * (prDistribution[index] - prIntervalMean);

        if (fabsf(diffRRarray[index] - diffMean) > (float) RRImaxThreshold) {
            RRImaxCount = RRImaxCount + 1;
        }

        index = index + 1;
    }//end


    float sumRR = 0;
    index = 0;
    while (index < (nbins - 1)) {
        sumRR = sumRR + rrArray[index];
        *RMSSD = *RMSSD + (rrArray[index + 1] - rrArray[index]) * (rrArray[index + 1] - rrArray[index]);
        index = index + 1;
    }
    float rrMean = sumRR / (float) nbins;
    *RMSSD = sqrtf(*RMSSD / (float) (nbins - 1));
    *RMSSD = *RMSSD / rrMean;

//   %%%%%%%%%%%%%  build hist %%%%%%%%%%%%%%%%%%%%%%
    float histPropability[histWidth] = {0};//zeros(1,histWidth);
    index = 0;
    //% coef=20;
    int num;
    float temp;
    int diffNum = length;//%%40;
    while (index < histNumber) {
        num = 0;
        while (num < diffNum) {
            temp = floorf(
                    (float) histNumber * (diffRRarray[num] - (float) deltaRRmin) / (float) (deltaRRmax - deltaRRmin));

            if (temp == (float) index) {
                histPropability[index] = histPropability[index] + 1;
            }//end
            num = num + 1;
        }// end
        index = index + 1;
    }//end
    float Shannon1 = 0;
    index = 0;
    float p;
    int histValidatedCount = 0;
    int histMaxCount = (int) histPropability[index];
    while (index < histWidth) {
        if (histPropability[index] == 0)
            p = 1;
        else {
            p = (float) histPropability[index] / (float) effectiveLength;//%histWidth;
            histValidatedCount = histValidatedCount + 1;
        }
        if (histPropability[index] > (float) histMaxCount) {
            histMaxCount = (int) histPropability[index];
        }//end

//            float tttt1=(1/(float)histNumber);
//            float tt2=log2((1/(float)histNumber));
//            float tt3=p*log2(p);
//            float tt4=tt3/tt2;

        Shannon1 = Shannon1 + (p * log2f(p)) / log2f(1 / (float) histNumber);
        index = index + 1;
    }// end

    // %  Shannon1=(-1)*Shannon1;
    *shannonEntropy = Shannon1;
}

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void findTruePeak(float *array, int rrValue, float *adjustedPeak) {

    int mFrom = rrValue - BEAT_MS10;
    if (mFrom < 0)
        mFrom = mFrom + ECG_BUFFER_LENGTH;
    // end

    int mTo = rrValue + BEAT_MS30;
    if (mTo > ECG_BUFFER_LENGTH) {
        mTo = mFrom - ECG_BUFFER_LENGTH;
    }//end

    int index = mFrom;
    int temPosition = rrValue;
    float tempMax = fabsf(array[rrValue]);

    while (index < mTo) {
        if ((fabsf(array[index])) > tempMax) {
            tempMax = (fabsf(array[index]));
            temPosition = index;
        }
        index = index + 1;
        if (index > ECG_BUFFER_LENGTH) {
            index = index - ECG_BUFFER_LENGTH;
        }
    }

    *adjustedPeak = (float) temPosition;

}

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int rWavePeak(float *beat, int rrValue, float Crit) {
    int positionTemp = rrValue;
    int sPosition = rrValue - BEAT_MS50;
    int tempLen = ECG_BUFFER_LENGTH;
    sPosition = counterAdjust(sPosition, tempLen);//%%%14:04 10/12 2018 ypz
    float sMax = fabsf(beat[sPosition]);
    int foundPosition = sPosition;
    int limit = counterAdjust(rrValue + BEAT_MS10 * 2, tempLen);//%%%

    while ((sPosition > 1) && (sPosition < limit)) {
        if (fabsf(beat[sPosition]) > sMax) {
            sMax = fabsf(beat[sPosition]);
            foundPosition = sPosition;
        }
        sPosition = sPosition + 1;
        if (sPosition == ECG_BUFFER_LENGTH) {
            sPosition = 0;
        }
    }

    if (foundPosition != rrValue) {
        rrValue = foundPosition;
    }

    // %  finish=rrValue;
    int start = rrValue;
    //  % start = counterAdjust(rrValue-MS120, ECG_BUFFER_LENGTH );%%%14:04 10/12 2018 ypz
    int finish = counterAdjust(rrValue - MS60, ECG_BUFFER_LENGTH);
    int maxPosition = rrValue;

    int t1, t2, t3, t4;

    while (start > finish)        // %%%14:04 10/12 2018 ypz
        //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    {
        t1 = counterAdjust(start - BEAT_MS6, ECG_BUFFER_LENGTH);
        t2 = counterAdjust(start + BEAT_MS6, ECG_BUFFER_LENGTH);
        t3 = counterAdjust(start - BEAT_MS10, ECG_BUFFER_LENGTH);
        t4 = counterAdjust(start + BEAT_MS10, ECG_BUFFER_LENGTH);

        if (((beat[start] - beat[t1]) > Crit) && ((beat[start] - beat[t2]) > Crit)) {
            if (((((beat[start] - beat[t3])) > 0) && ((beat[start] - beat[t4]) > 0)) ||
                (((beat[start] - beat[t3]) < 0) && (beat[start] - beat[t4] < 0))) {
                //  %  if(ECGBuffer(start)>0)
                maxPosition = start;
                break;
                // % end
            }//end
        }//end
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // % start=start+1;
        start = start - 1;
        if (start < 0)  //%%% 14:04 10/12 2018 ypz
        {
            start = start + ECG_BUFFER_LENGTH;
        }////end
    }// end

    if ((start > finish)) {
        int complete = counterAdjust(start + BEAT_MS10, ECG_BUFFER_LENGTH);
        int positionFind = start;
        float maxValue = beat[start];
        int index = counterAdjust(start - BEAT_MS10, ECG_BUFFER_LENGTH);
        while (index <= complete) {
            if (beat[index] > maxValue) {
                maxValue = beat[index];
                positionFind = index;
            }
            index = counterAdjust(index + 1, ECG_BUFFER_LENGTH);
        }
        if (positionFind != start) {
            start = positionFind;
        }

        if ((start != rrValue) && ((beat[maxPosition]) > fabsf(beat[positionTemp] / 12))) {
            rrValue = positionFind;//%maxPosition;
        }//end
    }//end

    return rrValue;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// int ecgMain::noiseFlagArray[3]={0};
// int ecgMain::noiseTag=0;
void QRSMorphology(QRSMorphology_CTX *ctx, float *beat, int flag, int *noiseFlag, int *morphologyType, int init) {

    // noiseFlag=0;
    // morphologyType=0;
    int *noiseFlagArray = ctx->noiseFlagArray;
    int *ctx_noiseTag = &(ctx->noiseTag);

    if (init) {
        (*ctx_noiseTag) = 0;
        memset(noiseFlagArray, 0, sizeof(int) * 3);
        return;
    }


    int index = 0;
    float QRSPeaks[10] = {0};
    int QRSPeaksPositions[10] = {0};

    int i = FIDMARK - BEAT_MS80;
    float maxQRS = beat[FIDMARK];
    float minQRS = beat[FIDMARK];

    while (i < FIDMARK + BEAT_MS160) {
        if ((beat[i] > maxQRS) && ((beat[i] > beat[i + 4]) && (beat[i] > beat[i - 4]))) {
            maxQRS = beat[i];
        } else {
            if (beat[i] < minQRS) {
                minQRS = beat[i];
            }
        }
        i = i + 1;
    }

    float Crit = 0.02f * (maxQRS - minQRS);
    float mainWaveHeight;

    if ((fabsf(maxQRS)) > fabsf(minQRS)) {
        mainWaveHeight = fabsf(maxQRS);
    } else {
        mainWaveHeight = fabsf(minQRS);
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // X=beat;
    float X[Fs];
    //       utility tempUtility;

    //  int  Fpa=200;
    //float returnValue=tempUtility.HFnoiseCheck(beat);
    float returnValue = HFnoiseCheck(beat);
    // tempUtility.LowPassFilter(beat,X);
    LowPassFilter(beat, X);
    int len = 90;
    float Xaux[90] = {0};
    int start = 80;
    int finish = 170;

    //tempUtility.ecgCopy(X,Xaux, start,finish);
    ecgCopy(X, Xaux, start, finish);
    //float mnoise=tempUtility.noiseLevel(X,len);
    float mnoise = noiseLevel(X, len);

//     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int leng = 90;
//        float xMax=tempUtility.ECGMax( Xaux,leng );
//        float xMin=tempUtility.ECGMin( Xaux,leng );
    // motionSum=sum(Xaux);
    float motionSum = ecgSum(Xaux, leng);
    //  float pp=motionSum/leng;
    start = 70;
    finish = 90;
    ecgCopy(X, Xaux, start, finish);
    // Xaux=X(70:160);
    float artifactsIndex = ecgSum(Xaux, leng);
    int kk = 0;
    while (kk < 2) {
        noiseFlagArray[kk] = noiseFlagArray[kk + 1];
        kk = kk + 1;
    }
    noiseFlagArray[2] = (*ctx_noiseTag);
    if (((noiseFlagArray[0] == 1) && (noiseFlagArray[1] == 1)) ||
        ((noiseFlagArray[1] == 1) && (noiseFlagArray[2] == 1)) ||
        ((noiseFlagArray[0] == 1) && (noiseFlagArray[2] == 1))) {
        if (mnoise > 9) {
            *noiseFlag = 1;
            mnoise = 33;
        }
    }

    if ((mnoise > 18) && (mnoise < 32) && (fabsf(minQRS) > maxQRS))
        mnoise = 0;
    //end
    if ((returnValue > 75) || (mnoise > 24) || ((artifactsIndex > 80) && (mnoise >
                                                                          19.1)))//%%||(artifactsIndex>110))%&&(mnoise>20)%%(mnoise>49))%%||indexNumber>=6||indexNumberRight>=6)     % Change from 61  46 to 51  35
        *noiseFlag = 1;
    //     %          SpectrumAnalysis( beat);
    if (mnoise > 18) {
        (*ctx_noiseTag) = 1;
        return;
    } else
        (*ctx_noiseTag) = 0;

    //          y=X;
    float temp;
    index = 0;
    int jj;
    i = (FIDMARK - BEAT_MS60);
    int peakIndex = 0;
    //      %   while ((i<QRSfinish+BEAT_MS50)&&((i+BEAT_MS10)<arrayLimit))
    while ((i < (BEAT_MS480))) {
        if ((i < (FIDMARK + 5)) && (i > 8)) {
            if ((fabsf(X[i] - X[i - BEAT_MS10]) > Crit / 2) && (fabsf(X[i] - X[i + BEAT_MS10]) > Crit / 2)) {
                if (((((X[i] - X[i - BEAT_MS10])) > 0) && ((X[i] - X[i + BEAT_MS10]) > 0)) ||
                    (((X[i] - X[i - BEAT_MS10]) < 0) && (X[i] - X[i + BEAT_MS10] < 0))) {
                    index = index + 1;

                    //        %%%%%%%%%%%%%%%%%%%%%%
                    kk = i;
                    if (X[kk] < (X[kk - BEAT_MS10] + X[kk + BEAT_MS10]) / 2) // %% downward
                    {
                        jj = kk - BEAT_MS10;
                        temp = X[jj];
                        peakIndex = jj;
                        jj = jj + 1;
                        while (jj < (kk + BEAT_MS10)) {
                            if (temp > X[jj]) {
                                temp = X[jj];
                                peakIndex = jj;
                            }//end
                            jj = jj + 1;
                        }//end
                    } else                         //  %%%%  upward
                    if ((X[kk] > (X[kk - BEAT_MS10] + X[kk + BEAT_MS10]) / 2) && (X[kk] > X[kk - 2 * BEAT_MS10])) {
                        jj = kk - BEAT_MS10;
                        temp = X[jj];
                        peakIndex = jj;
                        jj = jj + 1;
                        while (jj < (kk + BEAT_MS10)) {
                            if (temp < X[jj]) {
                                temp = X[jj];
                                peakIndex = jj;
                            }//end
                            jj = jj + 1;
                        }//end
                    }// end

                    //end

                    //         %%%%%%%%%%%%%%%%%% mainWaveHeight
                    //        %       if(abs(y(peakIndex))>mainWaveHeight*2/5)
                    if ((peakIndex >= 1) && (fabsf(X[peakIndex]) > mainWaveHeight * 1 / 12))//% 22:38 28/07 2018
                    {
                        QRSPeaksPositions[index - 1] = peakIndex;
                        QRSPeaks[index - 1] = temp;
                    } else {
                        index = index - 1;
                    }// end
                    //     %                i=peakIndex+2*BEAT_MS10;
                    i = i + 2 * BEAT_MS10;
                }// end
            }//end
        } else //%%%%%%%%%
        {
            //                %%%% i>FIDMARK
            //      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ((i > 4) && (i >= (FIDMARK - BEAT_MS40)) && (i <= (FIDMARK + BEAT_MS40))) {
                if ((fabsf(X[i] - X[i - BEAT_MS8]) > Crit / 2) && (fabsf(X[i] - X[i + BEAT_MS8]) > Crit / 2)) {
                    if (((((X[i] - X[i - BEAT_MS8])) > 0) && ((X[i] - X[i + BEAT_MS8]) > 0)) ||
                        (((X[i] - X[i - BEAT_MS8]) < 0) && (X[i] - X[i + BEAT_MS8] < 0))) {
                        index = index + 1;
                        //           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        kk = i;
                        if (X[kk] < (X[kk - BEAT_MS8] + X[kk + BEAT_MS8]) / 2) // %% downward
                        {
                            jj = kk - BEAT_MS8;
                            temp = X[jj];
                            peakIndex = jj;
                            jj = jj + 1;
                            while (jj < (kk + BEAT_MS8)) {
                                if (temp > X[jj]) {
                                    temp = X[jj];
                                    peakIndex = jj;
                                }//end
                                jj = jj + 1;
                            }//end
                        } else {//   %%%%  upward
                            if ((X[kk] > (X[kk - BEAT_MS8] + X[kk + BEAT_MS8]) / 2) &&
                                (X[kk] > X[kk + BEAT_MS20]))//  %% avoid little ripple
                            {
                                jj = kk - BEAT_MS8;
                                temp = X[jj];
                                peakIndex = jj;
                                jj = jj + 1;
                                while (jj < (kk + BEAT_MS8)) {
                                    if (temp < X[jj]) {
                                        temp = X[jj];
                                        peakIndex = jj;
                                    }//end
                                    jj = jj + 1;
                                }//end
                            }// end

                        }// end    else
                        //       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if ((index > 1) && (fabsf(X[peakIndex]) >= mainWaveHeight * 1 / 8) &&
                            (peakIndex != QRSPeaksPositions[index - 1]))//%0.58)
                        {
                            QRSPeaksPositions[index - 1] = peakIndex;
                            QRSPeaks[index - 1] = temp;
                        } else {
                            index = index - 1;
                        }//end
                        //     %               i=peakIndex+2*BEAT_MS6;
                        i = i + 2 * BEAT_MS6;
                    }//end    if(((((X[i]-X[i-BEAT_MS8]))>0)&&((X[i]-X[i+BEAT_MS8])>0))||(((X[i]-X[i-BEAT_MS8])<0)&&(X[i]-X[i+BEAT_MS8]<0)))
                }//end  if((fabs(X[i]-X[i-BEAT_MS8])>Crit/2)&&(fabs(X[i]-X[i+BEAT_MS8])>Crit/2))

                //     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            } else //%%%%%%%% end for  if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
            {
                //      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if ((i > BEAT_MS10) && (fabsf(X[i] - X[i - 2 * BEAT_MS10]) > Crit / 2) &&
                    (fabsf(X[i] - X[i + 2 * BEAT_MS10]) > Crit / 2)) {
                    if (((((X[i] - X[i - BEAT_MS10])) > 0) && ((X[i] - X[i + BEAT_MS10]) > 0)) ||
                        (((X[i] - X[i - BEAT_MS10]) < 0) && (X[i] - X[i + BEAT_MS10] < 0))) {
                        index = index + 1;
                        //        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        kk = i;
                        if ((X[kk] < (X[kk - 2 * BEAT_MS10] + X[kk + 2 * BEAT_MS10]) / 2) &&
                            (X[kk] < X[kk - 2 * BEAT_MS10])) //       %% downward
                        {
                            jj = kk - BEAT_MS10;
                            temp = X[jj];
                            peakIndex = jj;
                            jj = jj + 1;
                            while (jj < (kk + BEAT_MS10)) {
                                if (temp > X[jj]) {
                                    temp = X[jj];
                                    peakIndex = jj;
                                }//end
                                jj = jj + 1;
                            }// end
                        } else {         // %%%%  upward
                            if ((X[kk] > (X[kk - 2 * BEAT_MS10] + X[kk + 2 * BEAT_MS10]) / 2) &&
                                (X[kk] > X[kk + 2 * BEAT_MS10])) // %% avoid little ripple
                            {
                                jj = kk - BEAT_MS10;
                                temp = X[jj];
                                peakIndex = jj;
                                jj = jj + 1;
                                while (jj < (kk + BEAT_MS10)) {
                                    if (temp < X[jj]) {
                                        temp = X[jj];
                                        peakIndex = jj;
                                    }//end
                                    jj = jj + 1;
                                }//end
                            }//end
                        }// end   else

                        //         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        //              %  if(abs(y(peakIndex))>abs(y(FIDMARK)/2))
                        if ((index > 1) && (fabsf(X[peakIndex]) >= mainWaveHeight * 1 /
                                                                   4))//&&(peakIndex!=QRSPeaksPositions[index-1]))//%0.58)  22:38  28/07 2018
                        {
                            if ((index >= 2) && (QRSPeaksPositions[index - 2] != peakIndex)) {
                                QRSPeaksPositions[index - 1] = peakIndex;
                                QRSPeaks[index - 1] = temp;

                            } else {
                                index = index - 1;
                            }
                        } else {
                            index = index - 1;
                        }// end
                        //     %            i=peakIndex+2*BEAT_MS10;
                        i = i + 2 * BEAT_MS10;
                    }//end  line 368  if(((((X[i]-X[i-BEAT_MS10]))>0)&&((X[i]-X[i+BEAT_MS10])>0))||(((X[i]-X[i-BEAT_MS10])<0)&&(X[i]-X[i+BEAT_MS10]<0)))
                }// end    line 366     if((i>BEAT_MS10)&&(fabs(X[i]-X[i-2*BEAT_MS10])>Crit/2)&&(fabs(X[i]-X[i+2*BEAT_MS10])>Crit/2))


            }// end   if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
            //     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        }//end   if((i<(FIDMARK+5))&&(i>8))   else
        i = i + 1;
    }//end while



    //       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //       %%%  type  1:  QS
    //       %%  type 2
    switch (index) {

        case 1: {
            if ((QRSPeaks[0] < 0) && (fabsf(QRSPeaks[0]) > 0.2 * ampValue)) {
                *morphologyType = 1;  //%%  "QS"
            } else {
                *morphologyType = 2; //%%%  "R'"
            }

            if ((QRSPeaks[0] > 0) && (fabsf(QRSPeaks[0]) > 0.2 * ampValue)) {
                *morphologyType = 2; //%%%  "R'"
            } // end
        }
            break;

        case 2: {
            if ((QRSPeaks[0] < 0) && (QRSPeaks[1] < 0) && (fabsf(QRSPeaks[1]) > 0.2 * ampValue))
                *morphologyType = 1;// %%  "QS"
            else {
                if ((QRSPeaks[0] > 0) && ((fabsf(QRSPeaks[1]) > 5 * QRSPeaks[0])) &&
                    (fabsf(QRSPeaks[1]) > 0.2 * ampValue) && (QRSPeaks[1] < 0)) {
                    *morphologyType = 1; //%%  "QS"
                }//end

            }
            if (((QRSPeaks[0] < 0) && (QRSPeaks[1] > 0) && (fabsf(QRSPeaks[0]) < 2 * (QRSPeaks[1]))) ||
                ((QRSPeaks[0] > 0) && (QRSPeaks[1] < 0) && (fabsf(QRSPeaks[1]) < (QRSPeaks[0])))) {
                *morphologyType = 2; //%%
            }
        }

            break;

        case 3: {

            if ((fabsf(QRSPeaks[1]) > 0.2 * ampValue) && ((QRSPeaks[1] < 0) || (QRSPeaks[1] > 0)) &&
                (fabsf(QRSPeaks[1]) > QRSPeaks[0]) && (fabsf(QRSPeaks[1]) > QRSPeaks[2])) {
                if ((QRSPeaksPositions[2] - QRSPeaksPositions[0]) < BEAT_MS100)
                    *morphologyType = 1; //%%  "QS"
                break;
            }

            if (((QRSPeaks[0] < 0) && (QRSPeaks[1] > 0) && (fabsf(QRSPeaks[0]) < 2 * (QRSPeaks[1]))) ||
                ((QRSPeaks[0] > 0) && (QRSPeaks[1] < 0) && (fabsf(QRSPeaks[1]) < (QRSPeaks[0])))) {
                *morphologyType = 2;// %%
            }

            if (((QRSPeaks[1] < 0) && (QRSPeaks[2] > 0) && (fabsf(QRSPeaks[1]) < 3 * (QRSPeaks[2]))) ||
                ((QRSPeaks[1] > 0) && (QRSPeaks[2] < 0) && (fabsf(QRSPeaks[2]) < (QRSPeaks[1])))) {
                *morphologyType = 2; //%%
            }

            if ((QRSPeaks[0] > 0) && (QRSPeaks[2] > 0) && (fabsf(QRSPeaks[1]) < 2 * (QRSPeaks[2]))) {
                *morphologyType = 2; //%%
            }

            if ((QRSPeaks[0] > 0) && (QRSPeaks[2] > 0) && (QRSPeaks[1] < 0)) {
                *morphologyType = 2;// %%
            }
        }
            break;

        case 4:
            //                      morphologyType=2; %%

            // otherwise
            *morphologyType = 2;

        default:;


    }//end

}


int
findSinusPWaves(float *array, int rr, int Interval, int currentOnset, int tOffset, float thresh, float currentPsum) {

    //  int rrPeriod=rr;
    int jjEnd;
    int jjStart;
    int arrayStart;
    int arrayEnd;
    //   if((rrPeriod>BEAT_MS1500)&&(rrPeriod<2*BEAT_MS1000)){
    jjEnd = Interval;
    jjStart = Interval - rr;
    if (jjStart < 0) {
        jjStart = jjStart + peakLocationMod;
    }

    arrayStart = jjStart + tOffset - FIDMARK;
    if (arrayStart > peakLocationMod) {
        arrayStart = arrayStart - peakLocationMod;
    }

    arrayEnd = jjEnd - (FIDMARK - currentOnset);
    if (arrayEnd < 0) {
        arrayEnd = arrayEnd + peakLocationMod;
    }

    int len = arrayEnd - arrayStart;
    if (len < 0) {
        len = len + peakLocationMod;
    }
    //
    int arrayLength = len;
    float ampliferFactor = 200;
    float peakValue[5] = {0};
    int peakPosition[5] = {0};

    // memory allocation
    //   float *mask = (float*)malloc(arrayLength*sizeof(float));
    float mask[3 * BEAT_MS1000] = {0};
//      if(mask==NULL)
//      {
//          arrayLength=len;
//      }
//      memset(mask,0,(arrayLength*sizeof(float)));
    // float mask[arrayLength];//zeros(1,arrayLength);
    // float Xpc[arrayLength];
    //    float *Xpc = (float*)malloc(arrayLength*sizeof(float));
    float Xpc[3 * BEAT_MS1000] = {0};
//      if(Xpc==NULL)
//      {
    int jj = arrayStart;
    int kk = 0;
    float tempArray[2 * BEAT_MS1000] = {0};

    while (kk < len) {
        tempArray[kk] = array[jj];
        kk = kk + 1;
        jj = jj + 1;
        if (jj >= peakLocationMod) {
            jj = jj - peakLocationMod;
        }
    }
    //Fs=500;
    int Fpa = 25;
    // qrsFeaturesClass sinusTemp;
    LowPassFilterB(array, Fs, Fpa, Xpc);

    int peakNumber = 0;
    //utility tempNoise;
    float ruido = noiseLevel(Xpc, arrayLength);
    arrayCopy(Xpc, tempArray, arrayLength);
    int start;
    if ((thresh > 2 * ruido) && (thresh > BEAT_MS20)) {
        start = 0;
        while (start < (arrayLength - 1)) {
            mask[start] = tempArray[start + 1] - tempArray[start];
            start = start + 1;
        }

        start = 0;
        while (start <= (arrayLength - 1)) {
            if (mask[start] > 0)
                mask[start] = 1;
            else if (mask[start] < 0)
                mask[start] = -1;
            //end
            //end
            start = start + 1;
        }

        start = 0;
        //  float y[arrayLength];
        //   float * y = (float *)malloc(arrayLength);
        arrayCopy(mask, Xpc, arrayLength);
        // y=mask;
        while (start < (arrayLength - 1)) {
            mask[start] = Xpc[start + 1] - Xpc[start];
            start = start + 1;
        }

        int pOnset;
        float temp;
        float tt;
        int maxPeakPosition;
        int firstFlag;
        float mm;
        start = 0;
        int extraLogic;
        int lastOne = 0;


        while ((start < 3 * BEAT_MS1000) && (start < (arrayLength - 1))) {
            firstFlag = 0;
            mm = fabsf(tempArray[(start) + 1]);
            extraLogic = ((((start - lastOne) > BEAT_MS200) && (firstFlag == 0)) ||
                          (((start - lastOne) > 2 * BEAT_MS200) && (firstFlag == 1)));
            if (((mask[start] == -2) || (mask[start] == 2)) && (mm > 10) && (mm > thresh * 0.4) &&
                (mm < thresh * 1.5) && (extraLogic == 1)) {
                firstFlag = 1;
                lastOne = start + 1;
                peakValue[peakNumber] = tempArray[(start) + 1];
                maxPeakPosition = (start + 1);
                peakPosition[peakNumber] = start + 1;
                start = start + BEAT_MS20;
                peakNumber = peakNumber + 1;

                // TwaveFeaturesGet pTemp;
                float Pmax;
                int PmaxPosition;
                float Pmin;
                int PminPosition;

                temp = (float) peakPosition[peakNumber - 1] - BEAT_MS80;
                if (temp < 0)
                    temp = 0;
                //end
                pOnset = 0;
                tt = (float) peakPosition[peakNumber - 1] -
                     BEAT_MS20;// localMaxAlgorithm( float* array,int begin,int finish ,float*tmax,int*maxPosition,float*tmin,int*minPosition );
                if (temp < tt) {
                    localMaxAlgorithmT(tempArray, (int) temp, peakPosition[peakNumber] - BEAT_MS20, &Pmax,
                                       &PmaxPosition,
                                       &Pmin, &PminPosition);
                    // [ Pmax,PmaxPosition,Pmin,PminPosition ] = localMaxAlgorithm( array,temp,peakPosition(peakNumber)-BEAT_MS20 );
                    pOnset = PmaxPosition;
                }

                // %%%%% ---- Check noise level ----
                int pOffset = 0;
                temp = (float) peakPosition[peakNumber - 1] + BEAT_MS80;
                if (temp > (float) arrayLength)
                    temp = (float) arrayLength - 1;
                //end
                tt = (float) (peakPosition[peakNumber - 1] + BEAT_MS20);
                if (tt < temp) {
                    //  [ Pmax,PmaxPosition,Pmin,PminPosition ] = localMaxAlgorithm( array,peakPosition(peakNumber)+BEAT_MS20,temp );
                    localMaxAlgorithmT(tempArray, peakPosition[peakNumber - 1] + BEAT_MS20, (int) temp, &Pmax,
                                       &PmaxPosition,
                                       &Pmin, &PminPosition);
                    pOffset = PmaxPosition;
                    // %%%%% ----  p  wave  Validation ----
                    float slope1 =
                            fabsf(tempArray[maxPeakPosition] - tempArray[maxPeakPosition - BEAT_MS10]) / ampliferFactor;
                    float slope2 =
                            fabsf(tempArray[maxPeakPosition] - tempArray[maxPeakPosition + BEAT_MS10]) / ampliferFactor;
                    float slope3 = fabsf(tempArray[maxPeakPosition] - tempArray[maxPeakPosition - 3 * BEAT_MS10]) /
                                   ampliferFactor;
                    float slope4 = fabsf(tempArray[maxPeakPosition] - tempArray[maxPeakPosition + 3 * BEAT_MS10]) /
                                   ampliferFactor;
                    int slopeLogic =
                            (slope1 > 0.0005) && (slope2 > 0.0005) && (slope3 > 0.0005 * 4) && (slope4 > 0.0005 * 4);
                    int pValidation = (pOffset > pOnset) && ((pOffset - pOnset) < BEAT_MS150) &&
                                      ((maxPeakPosition - pOnset) > BEAT_MS30);
                    int pMagnitude = (fabsf(tempArray[maxPeakPosition] - tempArray[pOnset]) > 2 * ruido) &&
                                     (fabsf(tempArray[maxPeakPosition] - tempArray[pOffset]) > 2 * ruido);
                    if ((slopeLogic == 0) || (pValidation == 0) || pMagnitude == 0)
                        if (peakNumber >= 1)
                            peakNumber = peakNumber - 1;//%%%
                }  //if(tt<temp)
            }//end  if(((mask[start]==-2)||(mask[start]==2))&&(mm>10)&&(mm>thresh*0.5)&&(mm<thresh*1.5)&&(extraLogic==1))

            start = start + 1;
        } //end     while (start<(arrayLength-1))
    }// end     if((thresh>2*ruido)&&(thresh>BEAT_MS20))


    int ii = 0;
    int kk1, kk2;//,kk3;
    while (ii < arrayLength) {
        kk1 = (int) mask[ii];
        kk2 = (int) Xpc[ii];
        //         kk3=y[ii];
        ii++;
    }


    //free(mask);
    //free(Xpc);
    //   free(y);

    return peakNumber;
}

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//  int ecgMain:: apcString[vTachyLength]={0};
//  int ecgMain:: apcBigeminy=0;
//  int ecgMain:: apcTrigeminy;
//  int ecgMain:: apcBigString[vTachyLength+1]={0};
//  int ecgMain:: apcTriString[vTachyLength*2]={0};
//  int ecgMain:: longestApcTachycardiaCount=0;
//  int ecgMain:: apcTachycardiaSegment=0;
//  int ecgMain:: apcPointerCount=0;
//  int ecgMain:: apcBiSegmentCount=0;
//  int ecgMain:: apcBiLongest=0;
//  int ecgMain:: fastestApcTachycardia=0;
//  int ecgMain:: doubleApcCount=0;
//  int ecgMain:: apcdoubleArray[vTachyLength]={0};
//  int ecgMain:: apcrrSum;
//  int ecgMain:: apcrrString[vTachyLength]={0};
//  int ecgMain:: apcbiFlag=0;
//  int ecgMain:: apcflag=0;
//  int ecgMain:: apcCount=0;
//  int ecgMain:: apctriCount=0;
//  int ecgMain:: apctriFisrFlag=0;
//  int ecgMain:: apctriSecongFlag=0;
//  int ecgMain:: apcTriSegmentCount=0;
//  int ecgMain:: apcTriLongest=0;
//  float ecgMain:: normalTemplate[BEAT_MS200]={0};


void atrialPrematureCount(atrialPrematureCount_CTX *ctx, int apcFlag, int rr, int *doubleAvc, int *apcBigeminyCount,
                          int *apcTrigeminyCount,
                          int *atrialTachycardiaConut, int *longestApcTachycardia, int *fastestApcTachy, int init) {

    int *apcString = ctx->apcString; // static int apcString[vTachyLength];
    int *apcBigString = ctx->apcBigString; // static int apcBigString[vTachyLength + 1];
    int *apcTriString = ctx->apcTriString; // static int apcTriString[vTachyLength * 2];
    int *apcdoubleArray = ctx->apcdoubleArray; // static int apcdoubleArray[vTachyLength];
    int *apcrrString = ctx->apcrrString; // static int apcrrString[vTachyLength];
    float *normalTemplate = ctx->normalTemplate; // static float normalTemplate[BEAT_MS200];

    int *ctx_apcBigeminy = &(ctx->apcBigeminy);
    int *ctx_apcTrigeminy = &(ctx->apcTrigeminy);
    int *ctx_longestApcTachycardiaCount = &(ctx->longestApcTachycardiaCount);
    int *ctx_apcTachycardiaSegment = &(ctx->apcTachycardiaSegment);
    int *ctx_apcPointerCount = &(ctx->apcPointerCount);
    int *ctx_apcBiSegmentCount = &(ctx->apcBiSegmentCount);
    int *ctx_apcBiLongest = &(ctx->apcBiLongest);
    int *ctx_fastestApcTachycardia = &(ctx->fastestApcTachycardia);
    int *ctx_doubleApcCount = &(ctx->doubleApcCount);
    int *ctx_apcrrSum = &(ctx->apcrrSum);
    int *ctx_apcbiFlag = &(ctx->apcbiFlag);
    int *ctx_apcflag = &(ctx->apcflag);
    int *ctx_apcCount = &(ctx->apcCount);
    int *ctx_apctriCount = &(ctx->apctriCount);
    int *ctx_apctriFisrFlag = &(ctx->apctriFisrFlag);
    int *ctx_apctriSecongFlag = &(ctx->apctriSecongFlag);
    int *ctx_apcTriSegmentCount = &(ctx->apcTriSegmentCount);
    int *ctx_apcTriLongest = &(ctx->apcTriLongest);

    if (init) {

        //int  apcString[vTachyLength]={0};
        memset(apcString, 0, sizeof(int) * vTachyLength);
        (*ctx_apcBigeminy) = 0;
        (*ctx_apcTrigeminy) = 0;
        //int  apcBigString[vTachyLength+1]={0};
        memset(apcBigString, 0, sizeof(int) * (vTachyLength + 1));
        // int  apcTriString[vTachyLength*2]={0};
        memset(apcTriString, 0, sizeof(int) * vTachyLength * 2);
        (*ctx_longestApcTachycardiaCount) = 0;
        (*ctx_apcTachycardiaSegment) = 0;
        (*ctx_apcPointerCount) = 0;
        (*ctx_apcBiSegmentCount) = 0;
        (*ctx_apcBiLongest) = 0;
        (*ctx_fastestApcTachycardia) = 0;
        (*ctx_doubleApcCount) = 0;
        // int  apcdoubleArray[vTachyLength]={0};
        memset(apcdoubleArray, 0, sizeof(int) * vTachyLength);
        (*ctx_apcrrSum) = 0;
        // int  apcrrString[vTachyLength]={0};
        memset(apcrrString, 0, sizeof(int) * vTachyLength);
        (*ctx_apcbiFlag) = 0;
        (*ctx_apcflag) = 0;
        (*ctx_apcCount) = 0;
        (*ctx_apctriCount) = 0;
        (*ctx_apctriFisrFlag) = 0;
        (*ctx_apctriSecongFlag) = 0;
        (*ctx_apcTriSegmentCount) = 0;
        (*ctx_apcTriLongest) = 0;
        //float  normalTemplate[BEAT_MS200]={0};
        memset(normalTemplate, 0, sizeof(int) * (BEAT_MS200));

        return;
    }



    //  %%%%%%% atrial  tachycardia count %%%%%%
    float temp;
    int ii;
    if ((*ctx_apcPointerCount) == 0) {
        ii = 0;
        while (ii < (vTachyLength - 1)) {
            apcString[ii] = apcString[ii + 1];
            apcrrString[ii] = apcrrString[ii + 1];
            ii = ii + 1;
        }
        if (apcFlag == 1) {
            apcString[vTachyLength - 1] = 1;
            apcrrString[vTachyLength - 1] = rr;
        } else {
            apcString[vTachyLength - 1] = 0;
            apcrrString[vTachyLength - 1] = 0;
        }
    } else {
        if ((*ctx_apcPointerCount) > (*ctx_longestApcTachycardiaCount)) {
            (*ctx_longestApcTachycardiaCount) = (*ctx_apcPointerCount);
        }//end
        if (apcFlag == 1) {
            // %%%%%% long tachycardia %%%%%%
            (*ctx_apcflag) = 1;
            (*ctx_apcPointerCount) = (*ctx_apcPointerCount) + 1;
            (*ctx_apcrrSum) = (*ctx_apcrrSum) + rr;
            *doubleAvc = (*ctx_doubleApcCount);
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( doublePvcCount  ,pvcBiSegmentCount,pvcTriSegmentCount,pvcTachycardiaSegment,longestPvcTachycardiaCount,fastestPvcTachycardia);
            //    ventricularDataCopy( doubleAvc  , apcBigeminyCount,apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia,fastestApcTachy,(*ctx_doubleApcCount), (*ctx_apcBigeminy),apcTrigeminyCount,atrialTachycardiaConut,longestApcTachycardia,fastestApcTachy);
            //    return;
            //    void ecgMain::ventricularDataCopy( int doublePvcCount  ,int pvcBiSegmentCount,int pvcTriSegmentCount,int pvcTachycardiaSegment,int longestPvcTachycardiaCount,int fastestPvcTachycardia,int*doublePvc, int*vpcBigeminyCount,int*vpcTrigeminyCount,int*ventricularTachycardiaConut,int*longestPvcTachycardia,int*fastestPvcTachy)

            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        } else {
            temp = (float) 2 * (float) (*ctx_apcrrSum) / (float) (*ctx_apcPointerCount);
            float period = 60000 / temp;
            if (period > (float) (*ctx_fastestApcTachycardia)) {
                (*ctx_fastestApcTachycardia) = (int) period;
            }//end
            (*ctx_apcPointerCount) = 0;
            (*ctx_apcTachycardiaSegment) = (*ctx_apcTachycardiaSegment) + 1;
            (*ctx_apcflag) = 0;
            ii = 0;
            while (ii < 3) {
                apcString[ii] = 0;
                apcrrString[ii] = 0;
                ii = ii + 1;
            }
            clearBuffer(apcTriString, apcdoubleArray, apcBigString);
            // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( doublePvcCount  ,pvcBiSegmentCount,pvcTriSegmentCount,pvcTachycardiaSegment,longestPvcTachycardiaCount,fastestPvcTachycardia);
            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        }//end
    }// end

    if (((*ctx_apcflag) == 0) && (apcString[0] == 1) && (apcString[1] == 1) && (apcString[2] == 1)) {
        (*ctx_apcPointerCount) = 3;
        (*ctx_apcrrSum) = apcrrString[0] + apcrrString[1] + apcrrString[2];
        // [doublePvc, vpcBigeminyCount,vpcTrigeminyCount,ventricularTachycardiaConut,longestPvcTachycardia,fastestPvcTachy] = ventricularDataCopy( doublePvcCount  ,pvcBiSegmentCount,pvcTriSegmentCount,pvcTachycardiaSegment,longestPvcTachycardiaCount,fastestPvcTachycardia);
        ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                            (*ctx_apcTachycardiaSegment),
                            (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                            apcBigeminyCount,
                            apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
        return;
    }//end

    // %%%%%%% atrial tachycardia count  end %%%%%%;

    //   %%%%%%%  bigeminy counter  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

    if ((*ctx_apcBigeminy) == 0) {
        ii = 0;
        while (ii < 4) {
            apcBigString[ii] = apcBigString[ii + 1];
            ii = ii + 1;
        }//end
        if (apcFlag == 1) {
            apcBigString[3] = 1;
        } else {
            apcBigString[3] = 0;
        }//end
    } else {
        (*ctx_apcCount) = (*ctx_apcCount) + 1;
        if (((*ctx_apcCount) == 1) && (apcFlag == 0)) {
            (*ctx_apcbiFlag) = 1;
            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        }

        if ((apcFlag == 1) && ((*ctx_apcbiFlag) == 1) && ((*ctx_apcCount) == 2)) {
            (*ctx_apcCount) = 0;
            (*ctx_apcbiFlag) = 0;
            (*ctx_apcBigeminy) = (*ctx_apcBigeminy) + 1;
            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        } else {
            (*ctx_apcBiSegmentCount) = (*ctx_apcBiSegmentCount) + 1;
            if ((*ctx_apcBiLongest) < (*ctx_apcBigeminy)) {
                (*ctx_apcBiLongest) = (*ctx_apcBigeminy);
            }//end
            (*ctx_apcBigeminy) = 0;
            (*ctx_apcCount) = 0;

            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        }//end
    }// end

    if ((apcBigString[0] == 0) && (apcBigString[1] == 1) && (apcBigString[2] == 0) && (apcBigString[3] == 1)) {
        (*ctx_apcBigeminy) = (*ctx_apcBigeminy) + 1;
        ii = 0;
        while (ii < 4) {
            apcBigString[ii] = 0;
            ii = ii + 1;
        }
        (*ctx_apcbiFlag) = 0;
        (*ctx_apcCount) = 0;
        ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                            (*ctx_apcTachycardiaSegment),
                            (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                            apcBigeminyCount,
                            apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
        return;
    }//end

    //%%%%%%% bigeminy counter  end  %%%

    ii = 0;
    while (ii < 2) {
        apcdoubleArray[ii] = apcdoubleArray[ii + 1];
        ii = ii + 1;
    }
    if (apcFlag == 1)
        apcdoubleArray[2] = 1;
    else {
        apcdoubleArray[2] = 0;
    }//end

    //    %%%%%%%%%%% double apc %%%%%%%%%%%%
    if ((apcdoubleArray[0] == 1) && (apcdoubleArray[1] == 1)) {
        (*ctx_doubleApcCount) = (*ctx_doubleApcCount) + 1;
        ii = 0;
        while (ii < 3) {
            apcdoubleArray[ii] = 0;
            ii = ii + 1;
        }//end
        ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                            (*ctx_apcTachycardiaSegment),
                            (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                            apcBigeminyCount,
                            apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
        return;
    }//end

    //   %%%%%%%%%%%%%%% Trigeminy %%%%%%%%%%%
    if ((*ctx_apcTrigeminy) == 0) {
        ii = 0;
        while (ii < 5) {
            apcTriString[ii] = apcTriString[ii + 1];
            ii = ii + 1;
        }
        if (apcFlag == 1)
            apcTriString[5] = 1;
        else {
            apcTriString[5] = 0;
        }//end
    } else {
        (*ctx_apctriCount) = (*ctx_apctriCount) + 1;
        if (((*ctx_apctriCount) == 1) && (apcFlag == 0)) {
            (*ctx_apctriFisrFlag) = 1;
            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        } else {
            if (((*ctx_apctriCount) == 2) && (apcFlag == 0)) {
                (*ctx_apctriSecongFlag) = 1;
                ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                    (*ctx_apcTachycardiaSegment),
                                    (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                    apcBigeminyCount,
                                    apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
                return;
            }//end
        }// end

        if (((*ctx_apctriFisrFlag) == 1) && ((*ctx_apctriSecongFlag) == 1) && (apcFlag == 1) &&
            ((*ctx_apctriCount) == 3)) {
            (*ctx_apcTrigeminy) = (*ctx_apcTrigeminy) + 1;
            (*ctx_apctriFisrFlag) = 0;
            (*ctx_apctriSecongFlag) = 0;
            (*ctx_apctriCount) = 0;
        } else {
            (*ctx_apcTriSegmentCount) = (*ctx_apcTriSegmentCount) + 1;
            if ((*ctx_apcTriLongest) < (*ctx_apcTrigeminy)) {
                (*ctx_apcTriLongest) = (*ctx_apcTrigeminy);
            }//end
            (*ctx_apcTrigeminy) = 0;
            (*ctx_apctriCount) = 0;
            ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy),
                                (*ctx_apcTachycardiaSegment),
                                (*ctx_longestApcTachycardiaCount), (*ctx_fastestApcTachycardia), doubleAvc,
                                apcBigeminyCount,
                                apcTrigeminyCount, atrialTachycardiaConut, longestApcTachycardia, fastestApcTachy);
            return;
        }//end


    }// end

    if ((apcTriString[0] == 0) && (apcTriString[1] == 0) && (apcTriString[2] == 1) && (apcTriString[3] == 0) &&
        (apcTriString[4] == 0) && (apcTriString[5] == 1)) {
        (*ctx_apcTrigeminy) = (*ctx_apcTrigeminy) + 1;
        ii = 0;
        while (ii < 6) {
            apcTriString[ii] = 0;
            ii = ii + 1;
        }//end
        (*ctx_apctriCount) = 0;
    }//end
    ventricularDataCopy((*ctx_doubleApcCount), (*ctx_apcBigeminy), (*ctx_apcTrigeminy), (*ctx_apcTachycardiaSegment),
                        (*ctx_longestApcTachycardiaCount),
                        (*ctx_fastestApcTachycardia), doubleAvc, apcBigeminyCount, apcTrigeminyCount,
                        atrialTachycardiaConut,
                        longestApcTachycardia, fastestApcTachy);
    //  %%%%%%%%%%%%%%% Trigeminy  end  %%%%%%%%%%%

}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
float calculateCrit(const float *beat, int rrValue) {
    // %   ECG_BUFFER_LENGTH=3000;
    // BEAT_MS100 =50;
    int i = rrValue - BEAT_MS100;
    if (i < 0)
        i = i + ECG_BUFFER_LENGTH;
    //end

    float maxQRS = beat[i];
    float minQRS = beat[i];
    int limit = rrValue + BEAT_MS100;
    if (limit > ECG_BUFFER_LENGTH)
        limit = limit - ECG_BUFFER_LENGTH;
    //end
    while (i < limit) {
        if (beat[i] > maxQRS)
            maxQRS = beat[i];
        else if (beat[i] < minQRS)
            minQRS = beat[i];
        //end
        //end
        i = i + 1;

        if (i++ == ECG_BUFFER_LENGTH)
            i = 0;

    }//end

    // Crit=0.02*(maxQRS-minQRS);
    return 0.02f * (maxQRS - minQRS);

}


//  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// int ecgMain::ECGBufferIndex=0;
// float ecgMain::ECGBuffer[ECG_BUFFER_LENGTH]={0};
// float ecgMain::BufferV1[ECG_BUFFER_LENGTH]={0};
// float ecgMain::rWaveLocation[ECG_BUFFER_LENGTH]={0};
// float ecgMain::rWaveLocationV1[ECG_BUFFER_LENGTH]={0};

// float ecgMain::lastBeat[sampleRate]={0};
// float ecgMain::secondBeat[sampleRate]={0};

// float  ecgMain::prDistribution[histLength]={0};
// float  ecgMain::rrAFArray[histLength]={0};

// int ecgMain::lastRR=0;
// int ecgMain::lastRRV1=0;
// int ecgMain::peakNumber=0;
// int ecgMain::indexFlag=0;
// int ecgMain::learningNumber=0;
// int ecgMain::learningFlag=0;
// int ecgMain::complementFlag=0;
// int ecgMain::afEttopicCount=0;
// int ecgMain::rrAFCount=0;
// int ecgMain::pvcAFlag=0;
// int ecgMain::apcAFlag=0;
// int ecgMain::apcAFlagCount=0;
// int ecgMain::belowMs200=0;
// int ecgMain::avrQRSRR=0;
// int ecgMain::noiseCount=0;
// int ecgMain:: ppIntervalTotal=0;
// int ecgMain:: ppIntervalCount=0;
// float ecgMain::currentAverSpeak=0;
// float ecgMain:: maxQRSAerea=maxInteger;

// int ecgMain:: rbbbCount=0;;
// int ecgMain::lbbbCount=0;
// int ecgMain::japc=0;
// int ecgMain::apc=0;
// int ecgMain::wpwCount=0;

// int ecgMain::doublePvc=0;
// int ecgMain::vpcBigeminyCount=0;
// int ecgMain::vpcTrigeminyCount=0;
// int ecgMain::ventricularTachycardiaCount=0;
// int ecgMain::longestPvcTachycardia=0;
// int ecgMain::fastestPvcTachy=0;

// int  ecgMain::doubleAvc=0;
// int  ecgMain::apcBigeminyCount=0;
// int  ecgMain::atrialTachycardiaConut=0;
// int  ecgMain::longestApcTachycardia=0;
// int  ecgMain::apcTrigeminyCount=0;
// int  ecgMain::fastestApcTachy=0;

// int ecgMain::ventrEscapeCount=0;
// int ecgMain::junctEscapeCount=0;
// int ecgMain::templateFlag=0;

// int ecgMain:: sinusArrhythmiaType=0;
// int ecgMain:: sinusIrregularityType=0;
// int ecgMain:: sinusArrestType=0;
// int ecgMain:: sinusTachycardiaType=0;
// int ecgMain:: sickSinusSyndromeType=0;
// int ecgMain:: sinusBradycardiaType=0;

// int ecgMain:: twoDgreeBlockNumber=0;

// float ecgMain:: avrQRSWidth=maxInteger;
// float ecgMain:: currentTArea=maxInteger;
// float ecgMain:: currentArea=maxInteger;
// float ecgMain:: currentPeak=maxInteger;
// float ecgMain::currentPsum=maxInteger;
// float ecgMain:: tHeightAverage=maxInteger;
// float ecgMain:: currentRRSMD=0;

// float ecgMain::qrsRR[templateLength]={0};
// float ecgMain::qrsWidth[templateLength]={0};
// float ecgMain::qrsSpeak[templateLength]={0};
// float ecgMain::qrsArea[templateLength]={0};
// float ecgMain::tWaveArea[templateLength]={0};
// float ecgMain::Pheight[templateLength]={0};
// float ecgMain::Parea[templateLength]={0};
// float ecgMain::tNormalAverage[templateLength]={0};
// float ecgMain:: pwaveTemplate[2*BEAT_MS50]={0};
// float ecgMain::apcProcessing[MS60 ]={0};
// float ecgMain::AFapcProcessing[MS60 ]={300};

// PvcLibInfor ecgMain:: pvcLibInfo={0,{0},{{0},{0}}};

// int ecgMain:: stChangeCount=0;

// int ecgMain::pvcLength=0;
// int ecgMain::currentWidthMax=0;
// int ecgMain::currentQRS=0;
// //int ecgMain::apc=0;
// float ecgMain::pvcSum=0;
// float ecgMain::pvcStandard=0;
// float ecgMain::pvcWidthAverage=0;
// int  ecgMain::pvc=0;
// PvcNodeStru  ecgMain::pvcArray[templateLength]={{0,0}};
// float ecgMain::pvcStandardDerivation=0;
// float ecgMain::pvcStandardWidth=0;
// float ecgMain:: pvcWidthMaximum=0;
// float ecgMain::pvcAreaMaximum=0;
// int ecgMain:: pvcTemplateCount=0;
// int ecgMain::ST_IschaemicEpisode[20]={0};

// int ecgMain:: pvcWidthArray[6]={0};
// int ecgMain::pvcQrsWidthArray[6]={0};
// int ecgMain::effecitiveAFcount=0;

//  int  ecgMain::peakNumberV1=0;
//  /////////////////////////////////
//int  ecgMain::testCount=0;//////////////////////////

/*





	 */

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int peakAdjust(float beat[]) {

    int BEAT_MS16 = 8;
    // int BEAT_MS10=5;
    // int FIDMARK=200;

    int j = FIDMARK - 2 * BEAT_MS16;//%FIDMARK-8;
    float temp = beat[j];
    int Qindex = 0;
    while (j <= (FIDMARK + 2 * BEAT_MS10)) {
        if (beat[FIDMARK] > 0) {
            if (((beat[j + 1])) > temp) {
                temp = beat[j + 1];
                Qindex = j + 1;
            }//end
        } else {
            if ((fabsf(beat[j + 1])) > temp) {
                temp = fabsf(beat[j + 1]);
                Qindex = j + 1;
            }//end
        }//end
        j = j + 1;
    }//end

    //RPeak=Qindex;
    return Qindex;
}//end

//    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
void multiLeadsFeaturesGet(multiLeadsFeaturesGet_CTX *ctx, float ECGBuffer[][ECG_BUFFER_LENGTH], int rrValue, int *onset, int *offset, int *QPeak,
                           float *speakheight, float *rweight, int *beatBegin, int *beatEnd, float *amp, int *noiseFlag,
                           float *isoelectricLevel, int *notchFlag, float *IIArea, int temp) {

    //    int    BEAT_SAMPLE_RATE=500;
    int lenth = 4;
    int dataLength = sampleRate;//
    float tempBeat[BEAT_SAMPLE_RATE];

    int multiLeadArray[4] = {1, 5, 8, 10};//II, AVF,   V3,   V5
    float onArray[4];
    float offsetArray[4];
    float peakHeightArray[4];
    int noiseArray[4] = {0};
    memset(tempBeat, 0, sizeof(float) * sampleRate);
    memset(onArray, 0, sizeof(float) * 4);
    memset(offsetArray, 0, sizeof(float) * 4);
    memset(peakHeightArray, 0, sizeof(float) * 4);

    BeatFeaturesGet_CTX *ctx_beatFeaturesGetCtx = &(ctx->beatFeaturesGetCtx);

//%%   Computing II, AVF,   V3,   V5
    int pp;
    int rr = temp;
    int Interval = rrValue;
    int j = Interval - (SAMPLE_RATE / BEAT_SAMPLE_RATE) * FIDMARK;
    if (j < 0) {
        j = j + ECG_BUFFER_LENGTH;
    }//end
    int k;
    int start = j;
    int index = 0;
    while (index < 4) {
        k = 0;
        j = start;
        while ((k < (SAMPLE_RATE / BEAT_SAMPLE_RATE) * BEATLGTH))//%%&&(j<=ECG_BUFFER_LENGTH))
        {
            pp = multiLeadArray[index];
            tempBeat[k] = ECGBuffer[pp][j];
            j = j + 1;
            if (j == maxQueueLength - 1)
                j = 0;
            //end
            k = k + 1;
        }//end
        // %%%%%%     rWaveLocation(peakNumber)= rWavePeak( ECGBuffer,rrValue ,  Crit,ECG_BUFFER_LENGTH) ;
        //int onset=0;
        //int offset=0;
        //int QPeak=0;
        int SPeak = 0;
        float isoLevel = 0;
        *beatBegin = 0;
        *beatEnd = 0;
        *amp = 0;
        *rweight = 0;
        *isoelectricLevel = 0;
        float rrPVCLogic = 0;
        float Crit = 0;
        int interval = 0;
        Crit = calculateCrit(tempBeat, FIDMARK);
        float floatrrPeak = (float) peakAdjust(tempBeat);

        // %%%%%%%
        //  [ onset,offset,QPeak,speakheight,rweight,beatBegin,beatEnd,amp,noiseFlag ,isoelectricLevel,notchFlag] = BeatFeaturesGet(tempBeat,onset,offset,QPeak,SPeak,rr,beatBegin,beatEnd,amp ,Crit,0);
        BeatFeaturesGet(ctx_beatFeaturesGetCtx, tempBeat, &SPeak, interval, onset, offset, QPeak, speakheight, rweight, noiseFlag,
                        isoelectricLevel, notchFlag, 0);

        onArray[index] = (float) *onset;
        offsetArray[index] = (float) *offset;
        peakHeightArray[index] = *rweight;
        noiseArray[index] = *noiseFlag;
        index = index + 1;

    }//end while(index<4)

    // %% qrs wave onset
    // dataArray  = DataSort( onArray,lenth);
    DataSort(onArray, lenth);
    //dataArray=onArray;
    *onset = (int) floorf((onArray[1] + onArray[2]) / 2);
    // %% p wave offset
    DataSort(offsetArray, lenth);
    *offset = (int) floorf((offsetArray[1] + offsetArray[2]) / 2);
    //   %% p wave peak
    DataSort(peakHeightArray, lenth);
    *rweight = (peakHeightArray[1] + peakHeightArray[2]) / 2;


    //%% compute the II lead area




    float sum = 0;
    if (*onset > 0) {
        k = rrValue - (200 - *onset);
        int Limit = rrValue + (*offset - 200);
        int arrlength = abs(*offset - *onset);
        int jj = 0;
        while (jj < arrlength)//%ECGBuffer(2,k)
        {
            if (ECGBuffer[1][k] < 0)
                sum = sum + fabsf(ECGBuffer[1][k]);
            else
                sum = sum + (ECGBuffer[1][k]);
            k = k + 1;
            if (k == ECG_BUFFER_LENGTH)
                k = 0;

            jj = jj + 1;
        }//end
    }


    if ((noiseArray[0] == 1) && ((noiseArray[1] == 1) || (noiseArray[2] == 1) || (noiseArray[3] == 1))) {
        *notchFlag = 1;

    };
    *IIArea = sum;

}
//   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ECG_BUFFER_LENGTH/CoefficentConst

void MultiEcgCopy(float array[][ECG_BUFFER_LENGTH], float *ecgArray, int channelNo, int start, int finish) {
    //maxQueueLength=4000;
    //int len = abs(finish - start);
    // if(len<0)
    //   len=len+maxQueueLength;
    //end

//ecgArray=zeros(1,len);
    int jj = 0;
    int index = start;
    while (jj < 500) {

        ecgArray[jj] = array[channelNo][index];
        jj = jj + 1;
        index = index + 1;
        if (index > (maxQueueLength - 1))
            index = 0;
        //end
    }//end

}//end

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int ArrhythmiaAnalysis20230311(ArrhythmiaAnalysis20230311_CTX *ctx, const float *array, int ChannelNo, int *rrCount,
                               arrhythmia_type *result, int init,
                               int test_data_no) {
    float morphologyIndex;
    float recentPrInterval;
    int ST_MorphologyType;
    int ST_Change;
    float ST_ElavatorAmph;
    int type;
    float amplifierGainConst;
    float leadStandard;
    float fuzzyNormalResult;
    float areaIndex;
    float fuzzyPVCResult;
    float fuzzyApcResult;

    recentPrInterval = 0;
    morphologyIndex = 0;

    ST_MorphologyType = 0;
    ST_Change = 0;
    ST_ElavatorAmph = 0;
    type = 0;

    amplifierGainConst = 200;
    leadStandard = 0.2f;
    fuzzyNormalResult = 0;
    fuzzyApcResult = 0;
    areaIndex = 0;
    fuzzyPVCResult = 0;

    // ///////////
    float *tempBeat = ctx->tempBeat; // static float tempBeat[sampleRate];
    float *newArray = ctx->newArray; // static float newArray[BEAT_MS200];
    int *ctx_pvcNumber = &(ctx->pvcNumber);
    float *ectopicArray = ctx->ectopicArray; // static float ectopicArray[ectopicArrayLength];
    float *tempectopicArray = ctx->tempectopicArray; // static float tempectopicArray[ectopicArrayLength];
    float *rrAFArray = ctx->rrAFArray; // static float rrAFArray[histLength];
    int *ECGBufferIndex = &(ctx->ECGBufferIndex);
    float (*ECGBuffer)[ECG_BUFFER_LENGTH] = ctx->ECGBuffer; //static float ECGBuffer[ChannelMax][ECG_BUFFER_LENGTH];
    float *BufferV1 = ctx->BufferV1; // static float BufferV1[ECG_BUFFER_LENGTH];
    float *rWaveLocation = ctx->rWaveLocation; // static float rWaveLocation[ECG_BUFFER_LENGTH / CoefficentConst];
    float *rWaveLocationV1 = ctx->rWaveLocationV1; // static float rWaveLocationV1[ECG_BUFFER_LENGTH / CoefficentConst];
    float *lastBeat = ctx->lastBeat; // static float lastBeat[sampleRate];
    float *secondBeat = ctx->secondBeat; // static float secondBeat[sampleRate];
    float *prDistribution = ctx->prDistribution; // static float prDistribution[histLength];

    int *ctx_lastRR = &(ctx->lastRR);
    int *ctx_lastRRV1 = &(ctx->lastRRV1);
    int *ctx_peakNumber = &(ctx->peakNumber);
    int *ctx_indexFlag = &(ctx->indexFlag);
    int *ctx_learningNumber = &(ctx->learningNumber);
    int *ctx_learningFlag = &(ctx->learningFlag);
    int *ctx_complementFlag = &(ctx->complementFlag);
    int *ctx_afEttopicCount = &(ctx->afEttopicCount);
    int *ctx_rrAFCount = &(ctx->rrAFCount);
    int *ctx_pvcAFlag = &(ctx->pvcAFlag);
    int *ctx_apcAFlag = &(ctx->apcAFlag);
    int *ctx_apcAFlagCount = &(ctx->apcAFlagCount);
    int *ctx_belowMs200 = &(ctx->belowMs200);
    int *ctx_avrQRSRR = &(ctx->avrQRSRR);
    int *ctx_noiseCount = &(ctx->noiseCount);
    int *ctx_ppIntervalTotal = &(ctx->ppIntervalTotal);
    int *ctx_ppIntervalCount = &(ctx->ppIntervalCount);
    float *ctx_currentAverSpeak = &(ctx->currentAverSpeak);
    float *ctx_maxQRSAerea = &(ctx->maxQRSAerea);
    float (*ECGBufferMulti)[ECG_BUFFER_LENGTH /
                            CoefficentConst] = ctx->ECGBufferMulti; //static float ECGBufferMulti[ChannelMax][ECG_BUFFER_LENGTH / CoefficentConst];
    int *globalIndex = ctx->globalIndex; // static int globalIndex[ChannelMax];

    int *ctx_rbbbCount = &(ctx->rbbbCount);
    int *ctx_lbbbCount = &(ctx->lbbbCount);
    int *ctx_japc = &(ctx->japc);
    int *ctx_apc = &(ctx->apc);
    int *ctx_wpwCount = &(ctx->wpwCount);


    int *ctx_doublePvc = &(ctx->doublePvc);
    int *ctx_vpcBigeminyCount = &(ctx->vpcBigeminyCount);
    int *ctx_vpcTrigeminyCount = &(ctx->vpcTrigeminyCount);
    int *ctx_ventricularTachycardiaCount = &(ctx->ventricularTachycardiaCount);
    int *ctx_longestPvcTachycardia = &(ctx->longestPvcTachycardia);
    int *ctx_fastestPvcTachy = &(ctx->fastestPvcTachy);
    int *ctx_doubleAvc = &(ctx->doubleAvc);
    int *ctx_apcBigeminyCount = &(ctx->apcBigeminyCount);
    int *ctx_atrialTachycardiaConut = &(ctx->atrialTachycardiaConut);
    int *ctx_longestApcTachycardia = &(ctx->longestApcTachycardia);
    int *ctx_apcTrigeminyCount = &(ctx->apcTrigeminyCount);
    int *ctx_fastestApcTachy = &(ctx->fastestApcTachy);
    int *ctx_ventrEscapeCount = &(ctx->ventrEscapeCount);
    int *ctx_junctEscapeCount = &(ctx->junctEscapeCount);
    int *ctx_templateFlag = &(ctx->templateFlag);
    int *ctx_sinusArrhythmiaType = &(ctx->sinusArrhythmiaType);
    int *ctx_sinusIrregularityType = &(ctx->sinusIrregularityType);
    int *ctx_sinusArrestType = &(ctx->sinusArrestType);
    int *ctx_sinusTachycardiaType = &(ctx->sinusTachycardiaType);
    int *ctx_sickSinusSyndromeType = &(ctx->sickSinusSyndromeType);
    int *ctx_sinusBradycardiaType = &(ctx->sinusBradycardiaType);
    int *ctx_twoDgreeBlockNumber = &(ctx->twoDgreeBlockNumber);
    float *ctx_avrQRSWidth = &(ctx->avrQRSWidth);
    float *ctx_currentTArea = &(ctx->currentTArea);
    float *ctx_currentArea = &(ctx->currentArea);
    float *ctx_currentPeak = &(ctx->currentPeak);
    float *ctx_currentPsum = &(ctx->currentPsum);
    float *ctx_tHeightAverage = &(ctx->tHeightAverage);
    float *ctx_currentRRSMD = &(ctx->currentRRSMD);

    float *qrsRR = ctx->qrsRR; // static float qrsRR[templateLength];
    float *qrsWidth = ctx->qrsWidth; // static float qrsWidth[templateLength];
    float *qrsSpeak = ctx->qrsSpeak; // static float qrsSpeak[templateLength];
    float *qrsArea = ctx->qrsArea; // static float qrsArea[templateLength];
    float *tWaveArea = ctx->tWaveArea; // static float tWaveArea[templateLength];
    float *Pheight = ctx->Pheight; // static float Pheight[templateLength];
    float *Parea = ctx->Parea; // static float Parea[templateLength];
    float *tNormalAverage = ctx->tNormalAverage; // static float tNormalAverage[templateLength];
    float *pwaveTemplate = ctx->pwaveTemplate; // static float pwaveTemplate[2 * BEAT_MS50];
    float *apcProcessing = ctx->apcProcessing; // static float apcProcessing[MS60];
    float *AFapcProcessing = ctx->AFapcProcessing; // static float AFapcProcessing[MS60];

    PvcLibInfor *ctx_pvcLibInfo = &(ctx->pvcLibInfo);
    int *ctx_stChangeCount = &(ctx->stChangeCount);
    int *ctx_pvcLength = &(ctx->pvcLength);
    int *ctx_currentWidthMax = &(ctx->currentWidthMax);
    int *ctx_currentQRS = &(ctx->currentQRS);
    int *ctx_pvc = &(ctx->pvc);
    int *ctx_pvcTemplateCount = &(ctx->pvcTemplateCount);
    float *ctx_pvcSum = &(ctx->pvcSum);
    float *ctx_pvcStandard = &(ctx->pvcStandard);
    float *ctx_pvcWidthAverage = &(ctx->pvcWidthAverage);
    float *ctx_pvcStandardDerivation = &(ctx->pvcStandardDerivation);
    float *ctx_pvcStandardWidth = &(ctx->pvcStandardWidth);
    float *ctx_pvcWidthMaximum = &(ctx->pvcWidthMaximum);
    float *ctx_pvcAreaMaximum = &(ctx->pvcAreaMaximum);

    PvcNodeStru *pvcArray = ctx->pvcArray; // static PvcNodeStru pvcArray[templateLength];
    int *ST_IschaemicEpisode = ctx->ST_IschaemicEpisode; // static int ST_IschaemicEpisode[20];
    int *pvcWidthArray = ctx->pvcWidthArray; // static int pvcWidthArray[6];
    int *pvcQrsWidthArray = ctx->pvcQrsWidthArray; // static int pvcQrsWidthArray[6];

    ecgMainwdata_t *qrsFeatures = ctx->qrsFeatures; // static ecgMainwdata_t qrsFeatures[ECG_BUFFER_LENGTH / CoefficentConst];
    float *normalTemplate = ctx->normalTemplate; // static float normalTemplate[BEAT_MS200];

    int *rrInteral = ctx->rrInteral; // static int rrInteral[ECG_BUFFER_LENGTH];
    ecgMainwdata_t *MultiQrstFeatures = ctx->MultiQrstFeatures; // static ecgMainwdata_t MultiQrstFeatures[ECG_BUFFER_LENGTH];

    int *rrWaveLocation = ctx->rrWaveLocation; // static int rrWaveLocation[ECG_BUFFER_LENGTH / CoefficentConst];
    int *ctx_effecitiveAFcount = &(ctx->effecitiveAFcount);
    int *ctx_qrsFeaturesCountPreSec = &(ctx->qrsFeaturesCountPreSec);
    int *ctx_qrsFeaturesCountPre = &(ctx->qrsFeaturesCountPre);
    int *ctx_qrsFeaturesCount = &(ctx->qrsFeaturesCount);
    int *ctx_afEntryFlag = &(ctx->afEntryFlag);
    int *ctx_reEntryFlag_II = &(ctx->reEntryFlag_II);
    unsigned int *ctx_reEctopicFlag = &(ctx->reEctopicFlag);
    int *ctx_lastValue = &(ctx->lastValue);
    int *ctx_rrCountP = &(ctx->rrCountP);
    int *ctx_indexRr = &(ctx->indexRr);
    int *ctx_MultiPeakNumber = &(ctx->MultiPeakNumber);
    int *ctx_n = &(ctx->n);
    int *ctx_sIndex = &(ctx->sIndex);
    int *ctx_testCount = &(ctx->testCount);


    // context  =======================

    IntegerLowPassMulti_CTX *integerLowPassMultiCtx = &(ctx->integerLowPassMultiCtx);
    sinusArrhythmia_CTX *sinusArrhythmiaCtx = &(ctx->sinusArrhythmiaCtx);
    firstDegreeAVBlock_CTX *firstDegreeAvBlockCtx = &(ctx->firstDegreeAvBlockCtx);
    ventricularPrematureCount_CTX *ventricularPrematureCountCtx = &(ctx->ventricularPrematureCountCtx);
    atrialPrematureCount_CTX *atrialPrematureCountCtx = &(ctx->atrialPrematureCountCtx);
    QRSMorphology_CTX *qrsMorphologyCtx = &(ctx->qrsMorphologyCtx);
    IIRBandFilter_CTX *iirBandFilterCtx = &(ctx->iirBandFilterCtx);
    baseLineRemoval_CTX *baseLineRemovalCtx = &(ctx->baseLineRemovalCtx);
    IntegerLowPass_CTX *integerLowPassCtx = &(ctx->integerLowPassCtx);
    qrsDetectMultiLeads_CTX *qrsDetectMultiLeadsCtx = &(ctx->qrsDetectMultiLeadsCtx);
    baseLineRemovalMulti_CTX *baseLineRemovalMultiCtx = &(ctx->baseLineRemovalMultiCtx);
    NotchFilter_50hz_CTX *notchFilter50HzCtx = &(ctx->notchFilter50HzCtx);
    multiLeadsFeaturesGet_CTX *multiLeadsFeaturesGetCtx = &(ctx->multiLeadsFeaturesGetCtx);

    // ================================

    if (init == 1) {
        (*ctx_n) = 0;
        (*ctx_sIndex) = 0;
        int pp = 0;
        float p = 0;
        sinusArrhythmia(sinusArrhythmiaCtx, 0, &pp, &pp, &pp, &pp, &pp, &pp, 1);
        firstDegreeAVBlock(firstDegreeAvBlockCtx, 0, 1);
        ventricularPrematureCount(ventricularPrematureCountCtx, 0, 0, &pp, &pp, &pp, &pp, &pp, &pp, 1);
        atrialPrematureCount(atrialPrematureCountCtx, 0, 0, &pp, &pp, &pp, &pp, &pp, &pp, 1);
        QRSMorphology(qrsMorphologyCtx, &p, 0, &pp, &pp, 1);

        IIRBandFilter(iirBandFilterCtx, 20, 1);
        baseLineRemoval(baseLineRemovalCtx, 20, 1);
        IntegerLowPass(integerLowPassCtx, 20, 1);
        qrsDetectMultiLeads(qrsDetectMultiLeadsCtx, 0, 0, 1);
        baseLineRemovalMulti(baseLineRemovalMultiCtx, 0, 0, 1);
        IntegerLowPassMulti(integerLowPassMultiCtx, 0, 0, 1);
        NotchFilter_50hz(notchFilter50HzCtx, 20, 1, 1);
//
        (*ctx_qrsFeaturesCountPreSec) = 0;
        (*ctx_qrsFeaturesCountPre) = 0;
        (*ctx_qrsFeaturesCount) = 0;
        int ii = 0;
        while (ii < ECG_BUFFER_LENGTH / CoefficentConst) {

            qrsFeatures[ii].area = 0;
            qrsFeatures[ii].belowMs200 = 0;
            qrsFeatures[ii].morphologyType = 0;
            qrsFeatures[ii].noiseCount = 0;
            qrsFeatures[ii].notchFlag = 0;

            qrsFeatures[ii].offset = 0;
            qrsFeatures[ii].onset = 0;
            qrsFeatures[ii].Poffset = 0;
            qrsFeatures[ii].Ponset = 0;
            qrsFeatures[ii].pPeak = 0;

            qrsFeatures[ii].pPeakPositionFirst = 0;
            qrsFeatures[ii].pPeakPositionSecond = 0;
            qrsFeatures[ii].ppInterval = 0;
            qrsFeatures[ii].prDiff = 0;
            qrsFeatures[ii].prInterval = 0;
            qrsFeatures[ii].Psum = 0;
            qrsFeatures[ii].pType = 0;

            qrsFeatures[ii].Pwidth = 0;

            qrsFeatures[ii].qPeak = 0;
            qrsFeatures[ii].qrsHeight = 0;
            qrsFeatures[ii].qrsMainWave = 0;
            qrsFeatures[ii].qrsWidth = 0;
            qrsFeatures[ii].qtInterval = 0;
            qrsFeatures[ii].ration = 0;
            qrsFeatures[ii].rr = 0;
            qrsFeatures[ii].rrDiff = 0;
            qrsFeatures[ii].RRI = 0;
            qrsFeatures[ii].rrMax = 0;
            qrsFeatures[ii].shannon = 0;
            qrsFeatures[ii].speak = 0;
            qrsFeatures[ii].stChange = 0;
            qrsFeatures[ii].stType = 0;
            qrsFeatures[ii].tAera = 0;
            qrsFeatures[ii].tFirstPeak = 0;
            qrsFeatures[ii].tHeight = 0;
            qrsFeatures[ii].tOffset = 0;
            qrsFeatures[ii].Tonset = 0;
            qrsFeatures[ii].tSecondPeak = 0;
            qrsFeatures[ii].type = 0;
            //	qrsFeatures[ii].
            //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            MultiQrstFeatures[ii].rr = 0;
            MultiQrstFeatures[ii].onset = 0;
            MultiQrstFeatures[ii].offset = 0;
            MultiQrstFeatures[ii].qrsWidth = 0;
            MultiQrstFeatures[ii].Ponset = 0;
            MultiQrstFeatures[ii].Poffset = 0;
            MultiQrstFeatures[ii].Pwidth = 0;
            MultiQrstFeatures[ii].qrsHeight = 0;
            MultiQrstFeatures[ii].pPeak = 0;
            MultiQrstFeatures[ii].pType = 0;
            MultiQrstFeatures[ii].stChange = 0;
            MultiQrstFeatures[ii].stType = 0;
            MultiQrstFeatures[ii].speak = 0;
            MultiQrstFeatures[ii].noiseFlag = 0;
            MultiQrstFeatures[ii].tHeight = 0;
            MultiQrstFeatures[ii].Tonset = 0;
            MultiQrstFeatures[ii].tOffset = 0;
            MultiQrstFeatures[ii].qtInterval = 0;
            MultiQrstFeatures[ii].area = 0;
            MultiQrstFeatures[ii].prInterval = 0;
            MultiQrstFeatures[ii].ppInterval = 0;
            MultiQrstFeatures[ii].qrsMainWave = 0;
            MultiQrstFeatures[ii].tAera = 0;
            MultiQrstFeatures[ii].Psum = 0;
            MultiQrstFeatures[ii].tFirstPeak = 0;
            MultiQrstFeatures[ii].tSecondPeak = 0;

            MultiQrstFeatures[ii].rrDiff = 0;
            MultiQrstFeatures[ii].prDiff = 0;
            MultiQrstFeatures[ii].rrMax = 0;
            MultiQrstFeatures[ii].shannon = 0;
            MultiQrstFeatures[ii].ration = 0;
            MultiQrstFeatures[ii].noiseCount = 0;
            MultiQrstFeatures[ii].belowMs200 = 0;
            MultiQrstFeatures[ii].morphologyType = 0;


            MultiQrstFeatures[ii].RRI = 0;
            MultiQrstFeatures[ii].type = 0;
            MultiQrstFeatures[ii].pPeakPositionFirst = 0;
            MultiQrstFeatures[ii].pPeakPositionSecond = 0;
            MultiQrstFeatures[ii].vfFlag = 0;
            MultiQrstFeatures[ii].qpeak = 0;
            //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            rrWaveLocation[ii] = 0;

            ii++;
        }
//			//////

        ii = 0;
        while (ii < 10) {
            (*ctx_pvcLibInfo).Sdeep[ii] = 0;
            ii++;
        }

        ii = 0;
        while (ii < 8) {
            int jj = 0;
            while (jj < 100) {
                (*ctx_pvcLibInfo).informArray[ii][jj] = 0;
                jj++;
            }
            ii++;
        }
        (*ctx_pvcLibInfo).queueLength = 0;

//		////////////////
        ii = 0;
        while (ii < templateLength) {
            pvcArray[ii].qrsArea = 0;
            pvcArray[ii].qrsWidth = 0;
            ii++;
        }
//

        (*ctx_pvcNumber) = 0;
        ECGBufferIndex = 0;
        memset(ECGBuffer, 0, sizeof(float) * ChannelMax);//ChannelMax
        // ECGBuffer[ChannelMax][ECG_BUFFER_LENGTH]={0};
        memset(tempBeat, 0, sizeof(float) * sampleRate);
        memset(newArray, 0, sizeof(float) * BEAT_MS200);
        memset(ectopicArray, 0, ectopicArrayLength * sizeof(float));
        memset(tempectopicArray, 0, ectopicArrayLength * sizeof(float));

        memset(normalTemplate, 0, sizeof(float) * BEAT_MS200);
        //ECGBuffer[ECG_BUFFER_LENGTH]={0};
        //memset(ECGBuffer,0,sizeof(int)*ECG_BUFFER_LENGTH);
        //BufferV1[ECG_BUFFER_LENGTH]={0};
        //memset(BufferV1,0,sizeof(int)*ECG_BUFFER_LENGTH);
        //rWaveLocation[ECG_BUFFER_LENGTH]={0};
        memset(rWaveLocation, 0, sizeof(int) * ECG_BUFFER_LENGTH / CoefficentConst);
        // rWaveLocationV1[ECG_BUFFER_LENGTH]={0};
        memset(rWaveLocationV1, 0, sizeof(int) * ECG_BUFFER_LENGTH / CoefficentConst);

        //lastBeat[sampleRate]={0};
        memset(lastBeat, 0, sizeof(int) * sampleRate);
        // secondBeat[sampleRate]={0};
        memset(secondBeat, 0, sizeof(int) * sampleRate);
        //prDistribution[histLength]={0};
        memset(prDistribution, 0, sizeof(int) * histLength);
        //rrAFArray[histLength]={0};
        memset(rrAFArray, 0, sizeof(int) * histLength);
        memset(rrAFArray, 0, histLength * sizeof(float));

        // qrsDetect(0,ECGBuffer,1);
        // qrsDetectV1(0,ECGBuffer,1);

        (*ctx_lastRR) = 0;
        (*ctx_lastRRV1) = 0;
        (*ctx_peakNumber) = 0;
        (*ctx_indexFlag) = 0;
        (*ctx_learningNumber) = 0;
        (*ctx_learningFlag) = 0;
        (*ctx_complementFlag) = 0;
        (*ctx_afEttopicCount) = 0;
        (*ctx_rrAFCount) = 0;
        (*ctx_pvcAFlag) = 0;
        (*ctx_apcAFlag) = 0;
        (*ctx_apcAFlagCount) = 0;
        (*ctx_belowMs200) = 0;
        (*ctx_avrQRSRR) = 0;
        (*ctx_noiseCount) = 0;
        (*ctx_ppIntervalTotal) = 0;
        (*ctx_ppIntervalCount) = 0;
        (*ctx_currentAverSpeak) = 0;
        (*ctx_maxQRSAerea) = maxInteger;

        (*ctx_rbbbCount) = 0;
        (*ctx_lbbbCount) = 0;
        (*ctx_japc) = 0;
        (*ctx_apc) = 0;
        (*ctx_wpwCount) = 0;

        (*ctx_doublePvc) = 0;
        (*ctx_vpcBigeminyCount) = 0;
        (*ctx_vpcTrigeminyCount) = 0;
        (*ctx_ventricularTachycardiaCount) = 0;
        (*ctx_longestPvcTachycardia) = 0;
        (*ctx_fastestPvcTachy) = 0;

        (*ctx_doubleAvc) = 0;
        (*ctx_apcBigeminyCount) = 0;
        (*ctx_atrialTachycardiaConut) = 0;
        (*ctx_longestApcTachycardia) = 0;
        (*ctx_apcTrigeminyCount) = 0;
        (*ctx_fastestApcTachy) = 0;

        (*ctx_ventrEscapeCount) = 0;
        (*ctx_junctEscapeCount) = 0;
        (*ctx_templateFlag) = 0;

        (*ctx_sinusArrhythmiaType) = 0;
        (*ctx_sinusIrregularityType) = 0;
        (*ctx_sinusArrestType) = 0;
        (*ctx_sinusTachycardiaType) = 0;
        (*ctx_sickSinusSyndromeType) = 0;
        (*ctx_sinusBradycardiaType) = 0;

        (*ctx_twoDgreeBlockNumber) = 0;

        (*ctx_avrQRSWidth) = maxInteger;
        (*ctx_currentTArea) = maxInteger;
        (*ctx_currentArea) = maxInteger;
        (*ctx_currentPeak) = maxInteger;
        (*ctx_currentPsum) = maxInteger;
        (*ctx_tHeightAverage) = maxInteger;
        (*ctx_currentRRSMD) = 0;

        //qrsRR[templateLength]={0};
        memset(qrsRR, 0, sizeof(float) * templateLength);
        //qrsWidth[templateLength]={0};
        memset(qrsWidth, 0, sizeof(float) * templateLength);
        //qrsSpeak[templateLength]={0};
        memset(qrsSpeak, 0, sizeof(float) * templateLength);
        //qrsArea[templateLength]={0};
        memset(qrsSpeak, 0, sizeof(float) * templateLength);
        //tWaveArea[templateLength]={0};
        memset(tWaveArea, 0, sizeof(float) * templateLength);
        // Pheight[templateLength]={0};
        memset(Pheight, 0, sizeof(float) * templateLength);
        // Parea[templateLength]={0};
        memset(Parea, 0, sizeof(float) * templateLength);
        //tNormalAverage[templateLength]={0};
        memset(tNormalAverage, 0, sizeof(float) * templateLength);

        //pwaveTemplate[2*BEAT_MS50]={0};
        memset(pwaveTemplate, 0, sizeof(int) * 2 * BEAT_MS50);
        //apcProcessing[MS60 ]={0};
        memset(apcProcessing, 0, sizeof(int) * MS60);
        //AFapcProcessing[MS60 ]={300};
        memset(AFapcProcessing, 0, sizeof(int) * MS60);
        PvcLibInfor pvcLibInfo = {0, {0}, {{0}, {0}}};

        memset(rrInteral, 0, sizeof(int) * ECG_BUFFER_LENGTH);
        //rrInteral=zeros(1,ECG_BUFFER_LENGTH);

        (*ctx_stChangeCount) = 0;
        (*ctx_pvcLength) = 0;
        (*ctx_currentWidthMax) = 0;
        (*ctx_currentQRS) = 0;
        //int ecgMain::(*ctx_apc)=0;
        (*ctx_pvcSum) = 0;
        (*ctx_pvcStandard) = 0;
        (*ctx_pvcWidthAverage) = 0;
        (*ctx_pvc) = 0;
        PvcNodeStru pvcArray[templateLength] = {{0, 0}};
        (*ctx_pvcStandardDerivation) = 0;
        (*ctx_pvcStandardWidth) = 0;
        (*ctx_pvcWidthMaximum) = 0;
        (*ctx_pvcAreaMaximum) = 0;
        (*ctx_pvcTemplateCount) = 0;
        //ST_IschaemicEpisode[20]={0};
        memset(ST_IschaemicEpisode, 0, sizeof(int) * 20);
        //	pvcWidthArray[6]={0};
        memset(pvcWidthArray, 0, sizeof(int) * 6);
        //pvcQrsWidthArray[6]={0};
        memset(pvcQrsWidthArray, 0, sizeof(int) * 6);
        (*ctx_effecitiveAFcount) = 0;
        //	peakNumberV1=0;
        (*ctx_MultiPeakNumber) = 0;

        (*ctx_afEntryFlag) = 0;
        (*ctx_reEntryFlag_II) = 0;
        (*ctx_reEctopicFlag) = 0;
        (*ctx_lastValue) = 0;
        for (int i = 0; i < ChannelMax; ++i) {
            globalIndex[i] = -1;
        }
        (*ctx_rrCountP) = 0;
        (*ctx_indexRr) = 0;
        /////////////////////////////////
        (*ctx_testCount) = 0;
        init_arrhythmia_type(result);

        return 0;

    }  //  end if(init==1)


    for (int jj = 0; jj < ChannelNo; jj++) {
        float temp = array[jj];
        float ecgSample = IntegerLowPassMulti(integerLowPassMultiCtx, temp, jj, 0);
        ecgSample = baseLineRemovalMulti(baseLineRemovalMultiCtx, ecgSample, jj, 0);
        // for 50hz removal
        ecgSample = NotchFilter_50hz(notchFilter50HzCtx, ecgSample, jj, 0);

//        // test code --------------------
//        // ------------------------------
        globalIndex[jj] = globalIndex[jj] + 1;
        if (globalIndex[jj] >= (maxQueueLength - 1)) {
            globalIndex[jj] = 0;
        }
        ECGBufferMulti[jj][globalIndex[jj]] = ecgSample;

        //   int rrValue=qrsDetect2023(ecgSample,jj,&ECGBuffer,0);//,ECGBuffer#include "qrsDetect20230211.h"
        int rrValue = qrsDetectMultiLeads(qrsDetectMultiLeadsCtx, ecgSample, jj, 0);

        if (rrValue != 0) {
            if ((rrValue != 0) && (abs((*ctx_lastValue) - rrValue) > BEAT_MS200)) {
                //test(rrCount)=rrValue;
                rrWaveLocation[(*ctx_MultiPeakNumber)] = rrValue;
                (*ctx_MultiPeakNumber) = (*ctx_MultiPeakNumber) + 1;
                *rrCount = *rrCount + 1;
                if ((*ctx_MultiPeakNumber) > (ECG_BUFFER_LENGTH - 1)) {
                    (*ctx_MultiPeakNumber) = 0;
                }
                (*ctx_lastValue) = rrValue;
            }      //15/07 2023
            else {

                return 0;
            }
            // secondBeat=lastBeat;
            //lastBeat=tempBeat;

            arrayCopy(lastBeat, secondBeat, sampleRate);
            arrayCopy(tempBeat, lastBeat, sampleRate);
            float Crit;
            Crit = calculateCrit(ECGBuffer[1], rrValue);
///////////////////////////////////////////////////////////////////////////////////////////////////

// %%% %%%%  Compute  rr interval   %%%%%%%

            if (*rrCount >= 3) {

                temp = (float) (rrWaveLocation[(*ctx_MultiPeakNumber) - 1] -
                                rrWaveLocation[(*ctx_MultiPeakNumber) - 2]);
                if (temp < 1)
                    temp = temp + ECG_BUFFER_LENGTH;
                //end
                rrInteral[(*ctx_indexRr)] = (int) temp;//%%abs(test(mm)-test(mm-1));
                (*ctx_indexRr) = (*ctx_indexRr) + 1;
                if ((*ctx_indexRr) > (ECG_BUFFER_LENGTH - 1)) {
                    (*ctx_indexRr) = 0;
                }// end
                int peakNumberPre = (*ctx_MultiPeakNumber) - 1;
                if (peakNumberPre < 0)
                    peakNumberPre = peakNumberPre + ECG_BUFFER_LENGTH /
                                                    CoefficentConst;                  //  %%%%%%%%%%%%%%%%%%  14:04 10/12 2018 ypz%%%%%%%%%%%%%%%%%%  14:04 10/12 2018 ypz
                // end
                int Interval = rrWaveLocation[peakNumberPre];//%   peakLocationMod       // %%%%%%%%%%%%%%%%%%  14:04 10/12 2018 yp



                // %%%%%   %%%%%%  qrs wave   %%%%%%%
                //onset=0;offset=0;QPeak=0;speakheight=0;rweight=0;beatBegin=0;beatEnd=0;amp=0;noiseFlag=0 ;isoelectricLevel=0;notchFlag=0;
                int onset = 0;
                int SPeak = 0;
                int offset = 0;
                int QPeak = 0;
                float speakheight = 0;
                float rweight = 0;
                int noiseFlag = 0;
                float isoelectricLevel = 0;
                int notchFlag = 0;
                float IIArea;
                int beatBegin;
                int beatEnd;
                float amp;
                //%% Initialize the parameters
                // MultiLeadsParametersInitial( MultiQrstFeatures ,(*ctx_MultiPeakNumber));
                //[onset,offset,QPeak,speakheight,rweight,beatBegin,beatEnd,amp,noiseFlag ,isoelectricLevel,notchFlag,IIArea] = multiLeadsFeaturesGet(ECGBuffer,rrValue,temp);
                multiLeadsFeaturesGet(multiLeadsFeaturesGetCtx, ECGBufferMulti, rrValue, &onset, &offset, &QPeak, &speakheight, &rweight,
                                      &beatBegin,
                                      &beatEnd, &amp, &noiseFlag, &isoelectricLevel, &notchFlag, &IIArea, (int) temp);
                int rrPosition = rrValue;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr = rrInteral[(*ctx_indexRr) - 1];
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].onset = onset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].offset = offset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsWidth = offset - onset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsHeight = ECGBufferMulti[1][rrPosition];
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].qPeak = (float) QPeak;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].notchFlag = notchFlag;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].area = IIArea; // II
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].noiseFlag = noiseFlag;

                if (noiseFlag > 0) {
                    qrsFeatures[(*ctx_MultiPeakNumber)].prInterval = 0;
                    qrsFeatures[(*ctx_MultiPeakNumber)].noiseFlag = 1;
                    qrsFeatures[(*ctx_MultiPeakNumber)].tOffset = 0;
                    qrsFeatures[(*ctx_MultiPeakNumber)].Tonset = 0;
                    qrsFeatures[(*ctx_MultiPeakNumber)].ppInterval = 0;
                    qrsFeatures[(*ctx_MultiPeakNumber)].Ponset;
                    qrsFeatures[(*ctx_MultiPeakNumber)].Poffset = 0;
                    qrsFeatures[(*ctx_MultiPeakNumber)].pType = 0;
                    (*ctx_noiseCount) = (*ctx_noiseCount) + 1;
                    return 0;
                }

                // %%%%%%%%%%%%%%%%  For VT/VF Detect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  16:18  19/11 2022 noiseFlag
/*        float dispalyArray[1400];
        int  rr=rrInteral[(*ctx_indexRr)-1];
        int noiseTempCount=0;
        int j=rrValue-1400;
        j = counterAdjust( j, ECG_BUFFER_LENGTH );//%%%14:04 10/12 2018 ypz
        int  k=0;
        // %dispalyArray=zeros(1,1400);
        while(k<1400)
				{
          dispalyArray[k]=ECGBuffer[2][j] ;
          j=j+1;
          if(j==ECG_BUFFER_LENGTH)
             j=0;
          //end
          k=k+1;
        }//end

       int kk=0;
   while(kk<vfArrayLength)
	 {
       rrArray[kk]= rrArray[kk+1];
       qrsWidthArray[kk]= qrsWidthArray[kk+1];
       noiseStatistics[kk]=noiseStatistics[kk+1];
       if(noiseStatistics(kk)>=1)
         noiseTempCount=noiseTempCount+1;
       //end
       kk=kk+1;
		 }//end

    rrArray(vfArrayLength)=rrInteral[(*ctx_indexRr)-1];
    qrsWidthArray[vfArrayLength]=offset-onset;
    noiseStatistics[vfArrayLength]=noiseFlag;

      rrMean=mean(rrArray);
      rrStd=std(rrArray);

      qrsWidthMean=mean(qrsWidthArray);
      qrsWidthStd=std(qrsWidthArray);

        if((noiseTempCount<=1)&&((qrsWidthMean>95)||(qrsWidthMean==0)||(((rrStd/rrMean)>0.33)&&(qrsWidthMean>80))||(rrMean<250)||(rr>1000)))
          llll=0;
          [ A2,baseLinePDF,vt_vfFlag ] = vtSpectrumAnalysis(dispalyArray);
             if((A2>=0.45)&&(baseLinePDF<0.45))
               qrsFeatures((*ctx_peakNumber)).vfFlag=vt_vfFlag;%%  1:vt  2: vf
               if(vt_vfFlag==1)
                   (*ctx_pvc)=(*ctx_pvc)+1;
               else
                  vfCount=vfCount+1;
               end
               continue;
             end
        end

*/
                // %%%%%%%%%%%%%%% End for VF/VT  subroutine %%%%%%%%%%%%%%%%%%%%




                //%%%%%%%%% for AF detect  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 15:00 19/11 2022
                int rr = rrInteral[(*ctx_indexRr) - 1];
                int peakTempBufferOne;
                int peakTempBufferTwo;
                int peakTempBufferThree;
                int ectopicFlag = 0;
                int ectopicLogicFirst;
                int ectopicLogicSecond;
                int ectopicLogicThird;
                if ((*ctx_MultiPeakNumber) > 8) {
                    peakTempBufferOne = ((*ctx_MultiPeakNumber) - 1);
                    if (peakTempBufferOne < 0) {
                        peakTempBufferOne = peakTempBufferOne + maxQueueLength;
                    }//end
                    peakTempBufferTwo = ((*ctx_MultiPeakNumber) - 2);
                    if (peakTempBufferTwo < 0) {
                        peakTempBufferTwo = peakTempBufferTwo + maxQueueLength;
                    }//end
                    peakTempBufferThree = ((*ctx_MultiPeakNumber) - 3);
                    if (peakTempBufferThree < 1) {
                        peakTempBufferThree = peakTempBufferThree + maxQueueLength;
                    }//end
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].RRI =
                            (float) ((MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr) -
                                     (MultiQrstFeatures[peakTempBufferOne].rr));
                    ectopicLogicFirst = (((float) MultiQrstFeatures[peakTempBufferTwo].rr /
                                          (float) MultiQrstFeatures[peakTempBufferThree].rr) <= 0.83);
                    ectopicLogicSecond = (
                            ((float) MultiQrstFeatures[peakTempBufferOne].rr /
                             (float) MultiQrstFeatures[peakTempBufferTwo].rr) >= 1.3);
                    ectopicLogicThird = (
                            ((float) MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr /
                             (float) MultiQrstFeatures[peakTempBufferOne].rr) <= 0.9);
                    ectopicFlag = ectopicLogicFirst && ectopicLogicSecond && ectopicLogicThird;
                }//end

                (*ctx_rrAFCount) = (*ctx_rrAFCount) + 1;
                if ((ectopicFlag == 1) || ((*ctx_complementFlag) == 1)) {
                    if ((*ctx_rrAFCount) > 2) {
                        rrAFArray[(*ctx_rrAFCount) - 1] = 0;
                        rrAFArray[(*ctx_rrAFCount) - 2] = 0;
                        (*ctx_afEttopicCount) = (*ctx_afEttopicCount) + 1;
                    }//end
                }//end
                if ((MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsWidth) < 80)
                    rrAFArray[(*ctx_rrAFCount)] = (float) rr;
                //end
                if (rr >= BEAT_MS1500)
                    rrAFArray[(*ctx_rrAFCount)] = 0;
                //end

                if (((*ctx_rrAFCount) > 2) && ((*ctx_pvcAFlag) == 1) && (ectopicFlag == 0)) {
                    (*ctx_pvcAFlag) = 0;
                    rrAFArray[(*ctx_rrAFCount) - 1] = 0;
                    rrAFArray[(*ctx_rrAFCount) - 2] = 0;
                    //%  prDistribution((*ctx_rrAFCount)-1)=0;   %%%  14:34 09/05 2019
                }//end

                if ((*ctx_apcAFlag) == 1) {
                    (*ctx_apcAFlagCount) = (*ctx_apcAFlagCount) + 1;
                    (*ctx_apcAFlag) = 0;
                }//end

                peakTempBufferOne = ((*ctx_MultiPeakNumber) - 1);                          // %%%%%%%%%%%%%%%%%%%
                if (peakTempBufferOne < 1)
                    peakTempBufferOne = peakTempBufferOne + maxQueueLength;
                //end


                if (((*ctx_MultiPeakNumber) > 4) && (MultiQrstFeatures[peakTempBufferOne].pPeakPositionFirst != 0)) {
                    prDistribution[(*ctx_rrAFCount)] =
                            FIDMARK - MultiQrstFeatures[peakTempBufferOne].pPeakPositionFirst;
                }
                //end

                if (rr < BEAT_MS200) {
                    (*ctx_belowMs200) = (*ctx_belowMs200) + 1;
                }
                //end
                if ((noiseFlag == 1) && ((*ctx_rrAFCount) > 1)) {
                    rrAFArray[(*ctx_rrAFCount) - 1] = 0;
                    rrAFArray[(*ctx_rrAFCount)] = 0;
                }// end

                if ((*ctx_rrAFCount) == histLength) {
                    temp = (float) (rr - (*ctx_avrQRSRR)) / (float) (*ctx_avrQRSRR);
                    if ((temp > -0.6) && (temp < -0.2)) {
                        rrAFArray[(*ctx_rrAFCount)] = 0;
                    }//end
                }//end

                float RMSSD = 0;
                float prIntervalDiff = 0;
                float RRImaxCount = 0;
                float shannonEntropy;

                if ((*ctx_rrAFCount) == histLength) {
                    shannonEntropyCompute(rrAFArray, prDistribution, &shannonEntropy, &RMSSD, &prIntervalDiff,
                                          &RRImaxCount);
                    (*ctx_currentRRSMD) = RMSSD;
/*             afSegmentCount=afSegmentCount+1;
             afResultArray(afSegmentCount).rrDiff=RMSSD;
             afResultArray(afSegmentCount).prDiff=prIntervalDiff;
             afResultArray(afSegmentCount).rrMax=RRImaxCount;
            afResultArray(afSegmentCount).shannon=shannonEntropy;
 */
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].rrDiff = (int) RMSSD;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].prDiff = (int) prIntervalDiff;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].shannon = shannonEntropy;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].rrMax = (int) RRImaxCount;
                    if ((RMSSD > 0.10) &&
                        ((shannonEntropy > 0.35) || ((prIntervalDiff > 25) && (shannonEntropy > 0.30)))) {
                        (*ctx_effecitiveAFcount) = (*ctx_effecitiveAFcount) + 1;
                        if ((*ctx_apc) < 0)
                            (*ctx_apc) = 0;
                        //end
                    }//end

                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].belowMs200 = (*ctx_belowMs200);

//%             if((shannonEntropy>0.3)&&(RMSSD>0.1))
//%             end
                    (*ctx_noiseCount) = 0;
                    (*ctx_rrAFCount) = 0;
                    (*ctx_belowMs200) = 0;
                    (*ctx_apcAFlagCount) = 0;
                    // MotionArtifactAreaCount=0;
                    // MotionArtifactArea=0;
                }//end

                //%%%%%% end AF Segment %%%%%%




                // %%%%%%%%%%% P Wave detect    %%%%%%%%%%%%%%%  pay more to preToffset
                {
                    int j = rrValue - (SAMPLE_RATE / BEAT_SAMPLE_RATE) * FIDMARK;
                    if (j < 0) {
                        j = j + ECG_BUFFER_LENGTH;
                    }
                    for (int i = 0; i < (SAMPLE_RATE / BEAT_SAMPLE_RATE) * BEATLGTH; i++) {
                        tempBeat[i] = ECGBufferMulti[1][(i + j) % ECG_BUFFER_LENGTH]; // get II lead value
                    }
                }

                //Ponset=0;pPeak=0;Poffset=0 ; pPeakPositionFirst=0;pPeakPositionSecond=0;type=0 ;Psum=0; preToffset=0;
                //[Ponset,pPeak,Poffset , pPeakPositionFirst,pPeakPositionSecond,type ,Psum] = multiLeadsPWaveFeaturesGet(ECGBuffer,onset,rrValue,isoelectricLevel,preToffset );
                //   $$$$$$$$$$$  p wave detect  $$$$$$$$$$
                int Ponset = 0;
                float pPeak = 0;
                int Poffset = 0;
                int pPeakPositionFirst = 0;
                int pPeakPositionSecond = 0;
                int type = 0;
                float Psum = 0.0f;
                //                 pWaveDetected pClassTemplate;
                PWaveletDetect(tempBeat, onset, rr, isoelectricLevel, &Ponset, &pPeak, &Poffset, &pPeakPositionFirst,
                               &pPeakPositionSecond, &type, &Psum, 0);


                int pvcFlag = 0;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].prInterval = 0;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].ppInterval = 0;
//%          qrsFeatures((*ctx_peakNumber)).rrMax=prShannon;   %%% ypz 11:15 10/05 2019
// %       qrsFeatures((*ctx_peakNumber)).tHeight=tHeight;
//%         qrsFeatures((*ctx_peakNumber)).tOffset=tOffset;
//%         qrsFeatures((*ctx_peakNumber)).tFirstPeak=tFirstPeak;
//%         qrsFeatures((*ctx_peakNumber)).tSecondPeak=tSecondPeak;
//%         qrsFeatures((*ctx_peakNumber)).tAera=tSum;
//%
//%          MultiQrstFeatures((*ctx_peakNumber)).rr=rr;

                MultiQrstFeatures[(*ctx_MultiPeakNumber)].Ponset = Ponset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].Poffset = Poffset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].Pwidth = Poffset - Ponset;
                if ((Ponset == 0) && (Poffset == 0)) {
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].pPeak = 0;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].pType = 0;
                } else {
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].pPeak = pPeak;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].pType = type;
                }//end
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].speak = fabsf(speakheight);

                MultiQrstFeatures[(*ctx_MultiPeakNumber)].pPeakPositionFirst = pPeakPositionFirst;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].pPeakPositionSecond = pPeakPositionSecond;//%tbSum
                // %     qrsFeatures((*ctx_peakNumber)).tbSum=tbSum;
                //%     qrsFeatures((*ctx_peakNumber)).offset=offset;
                MultiQrstFeatures[(*ctx_MultiPeakNumber)].onset = onset;
                if ((Ponset > 1) && (onset > 1)) {
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].prInterval = onset - Ponset;
                } else {
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].prInterval = 0;
                }//end
                int peakTempBuffer;
                if (type != 0) {
                    peakTempBuffer = ((*ctx_MultiPeakNumber) - 1);
                    if (peakTempBuffer < 1) {
                        peakTempBuffer = peakTempBuffer + maxQueueLength;
                    }//end
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].ppInterval =
                            rrInteral[(*ctx_indexRr)] + Ponset - MultiQrstFeatures[peakTempBuffer].Ponset;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].prInterval = onset - Ponset;
                    (*ctx_ppIntervalTotal) =
                            (*ctx_ppIntervalTotal) + rrInteral[(*ctx_indexRr)] + Ponset -
                            MultiQrstFeatures[peakTempBuffer].Ponset;
                    (*ctx_ppIntervalCount) = (*ctx_ppIntervalCount) + 1;
                } else {
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].ppInterval = 0;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].prInterval = 0;
                }//end

                MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsMainWave = rweight;

                MultiQrstFeatures[(*ctx_MultiPeakNumber)].Psum = Psum;



// %%%%%% T Wave detect  %%%%%%%%%%%%%%%%%
                if (*rrCount > 3) {
                    // %  if(qrsFeatures(peakNumberPre).noiseFlag==0)

                    twav_data_t tPoint;
                    tPoint.tFirstPeak = 0;
                    tPoint.tSecondPeak = 0;
                    tPoint.type = 0;
                    tPoint.tHeight = 0;
                    tPoint.Tonset = 0;
                    tPoint.tOffset = 0;
                    tPoint.tAera = 0;
                    tPoint.qtInterval = 0;

                    // tFirstPeak=0;tSecondPeak=0;type=0;tHeight=0;Tonset=0;tOffset=0;tSum=0;
                    // lastPeak=(*ctx_MultiPeakNumber)-1;
                    // if(lastPeak<1)
                    //     lastPeak=lastPeak+ECG_BUFFER_LENGTH;
                    // end
                    float SPeak = 0;
                    float isoelectricLevel = 0;
                    float tLimit = 0;
                    float Crit = 0;
                    // [tFirstPeak,tSecondPeak,type,tHeight,Tonset,tOffset,tSum] = multiLeadsTWaveFeaturesGet(ECGBuffer, MultiQrstFeatures(lastPeak).onset, (*ctx_lastValue), MultiQrstFeatures(lastPeak).offset, SPeak,isoelectricLevel ,tLimit,Crit);
                    // TwaveFeaturesGet tClassTemplate;
                    int peakNumberPre = (*ctx_MultiPeakNumber) - 1;
                    if (peakNumberPre < 0) {
                        peakNumberPre = peakNumberPre + ECG_BUFFER_LENGTH;
                    }
                    tWavesFeaturesGet(lastBeat, MultiQrstFeatures[peakNumberPre].onset, rr,
                                      MultiQrstFeatures[peakNumberPre].offset,
                                      (int) MultiQrstFeatures[peakNumberPre].speak, isoelectricLevel, tLimit, Crit,
                                      &tPoint);

                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].tHeight = tPoint.tHeight;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].tOffset = tPoint.tOffset;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].qtInterval = tPoint.tOffset;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].tFirstPeak = tPoint.tFirstPeak;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].tSecondPeak = tPoint.tSecondPeak;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].tAera = tPoint.tAera;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].type = type;
                    MultiQrstFeatures[(*ctx_MultiPeakNumber)].Tonset = tPoint.Tonset;

                }//end   T Wave detect


                //%%%%%%%   T Wave detect   end  %%%%%%%%%%%%%%%%%%%%%%%%%%%

                //     %    (*ctx_MultiPeakNumber)=(*ctx_MultiPeakNumber)+1;
//%%%%%% 26/11  2022     begin

                peakNumberPre = (*ctx_MultiPeakNumber) - 1;
                if (peakNumberPre < 0) {
                    peakNumberPre = peakNumberPre +
                                    peakLocationMod;                             //  %%%%%%%%%%%%%%%%%%  14:04 10/12 2018 ypz
                }//en
                if (((*ctx_learningNumber) < templateLength) && ((*ctx_MultiPeakNumber) > 4)) {

                    qrsRR[(*ctx_learningNumber)] = (float) rrInteral[(*ctx_indexRr) - 1];
                    qrsWidth[(*ctx_learningNumber)] = (float) (offset - onset);
                    // %peakNumberPre=(*ctx_MultiPeakNumber)-1;
                    qrsSpeak[(*ctx_learningNumber)] = fabsf(MultiQrstFeatures[peakNumberPre].speak);
                    qrsArea[(*ctx_learningNumber)] = fabsf(MultiQrstFeatures[peakNumberPre].area);
                    tWaveArea[(*ctx_learningNumber)] = (MultiQrstFeatures[peakNumberPre].tAera);
                    Pheight[(*ctx_learningNumber)] = (MultiQrstFeatures[peakNumberPre].pPeak);//%
                    Parea[(*ctx_learningNumber)] = (MultiQrstFeatures[peakNumberPre].Psum);
                    tNormalAverage[(*ctx_learningNumber)] = fabsf(MultiQrstFeatures[peakNumberPre].tHeight);
                    (*ctx_learningNumber) = (*ctx_learningNumber) + 1;
                    if ((*ctx_learningNumber) >= templateLength) {
                        float AverResult = templateCompute(qrsRR, templateLength);
                        (*ctx_avrQRSRR) = (int) AverResult;
                        float currentRR = (float) (*ctx_avrQRSRR);
                        AverResult = templateCompute(qrsWidth, templateLength);
                        (*ctx_avrQRSWidth) = AverResult;
                        (*ctx_currentTArea) = templateCompute(tWaveArea, templateLength);
                        (*ctx_currentAverSpeak) = templateCompute(qrsSpeak, templateLength);
                        (*ctx_currentArea) = templateCompute(qrsArea, templateLength);
                        (*ctx_currentPeak) = templateCompute(Pheight, templateLength);
                        (*ctx_currentPsum) = templateCompute(Parea, templateLength);
                        (*ctx_tHeightAverage) = templateCompute(tNormalAverage, templateLength);
                        (*ctx_learningFlag) = 1;
                    }//end
                }//end
// %%%%%%%%%%%%
                float qrsWidthNarrow = 0;
                float qrsWidthWide = 0;
                float qrsIntervalStrong = 0;
                float qrsIntervalWeak = 0;
                float qrsComplementWeak = 0;
                float qrsComplementStrong = 0;
                float inputData;
                if ((peakNumberPre >= 1) && (((MultiQrstFeatures[peakNumberPre].qrsWidth)) >= 10) &&
                    ((MultiQrstFeatures[peakNumberPre].qrsWidth) <= 300)) {
                    inputData = (float) (MultiQrstFeatures[peakNumberPre].qrsWidth * 2);

                    if ((inputData < 268) && (inputData > 200)) inputData = 200;
                    qrsWidthNarrow = (float) FuzzyReasoning(1, inputData);
                    qrsWidthWide = (float) FuzzyReasoning(2, inputData);

                }//end

                if ((*ctx_complementFlag) == 3) (*ctx_complementFlag) = 0;

                if ((*ctx_complementFlag) >= 1) (*ctx_complementFlag) = (*ctx_complementFlag) + 1;

                // %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15:35 26/11  2022
                int peakLoc;
                if ((*ctx_MultiPeakNumber) >= 8) {
                    peakLoc = (*ctx_MultiPeakNumber) - 2;
                    if (peakLoc < 0) {
                        peakLoc = peakLoc + peakLocationMod;
                    }//end
                    int inputTemp = 0;
                    if (((*ctx_MultiPeakNumber) >= 3) && ((*ctx_complementFlag) == 0) &&
                        ((MultiQrstFeatures[peakLoc].rr) > (*ctx_avrQRSRR))) {
                        if (((float) MultiQrstFeatures[peakLoc].rr / (float) (*ctx_avrQRSRR)) > 1.2f) {  //for huge rr
                            inputData = (float) (MultiQrstFeatures[peakNumberPre].rr - (*ctx_avrQRSRR)) /
                                        (float) (*ctx_avrQRSRR);
                        } else {
                            peakTempBuffer = ((*ctx_MultiPeakNumber) - 2);
                            if (peakTempBuffer < 1) peakTempBuffer = peakTempBuffer + maxQueueLength;

                            inputData = (float) (MultiQrstFeatures[peakNumberPre].rr - MultiQrstFeatures[peakLoc].rr) /
                                        (float) MultiQrstFeatures[peakTempBuffer].rr;
                        }//end
                    } else {
                        //%%%%%%%%%%%%%% Continuous PVC   // %%%%%%%%%%%%%%%%%%
                        peakTempBuffer = ((*ctx_MultiPeakNumber) -
                                          1);                              // %%%%%%%%%%%%%%%%%%%%%
                        if (peakTempBuffer < 0) {
                            peakTempBuffer = peakTempBuffer + maxQueueLength;
                        }//end
                        float rrTempInterval = (float) (MultiQrstFeatures[peakLoc].rr) /
                                               (float) (MultiQrstFeatures[peakTempBuffer].rr);
                        // %   if(((*ctx_complementFlag)==2)&&((qrsFeatures((*ctx_peakNumber)-2).rr)/(qrsFeatures((*ctx_peakNumber)-1).rr)>0.9)&&((qrsFeatures((*ctx_peakNumber)-2).rr)/(qrsFeatures((*ctx_peakNumber)-1).rr)<1.17))
                        if (((((*ctx_complementFlag) == 2) && ((rrTempInterval > 0.9) && (rrTempInterval < 1.17))) ||
                             rrTempInterval < 0.88) || rrTempInterval > 1.2) {
                            inputData = (float) (MultiQrstFeatures[peakNumberPre].rr - (*ctx_avrQRSRR)) /
                                        (float) (*ctx_avrQRSRR);//%%%%% continuous PVC (*ctx_apc)
                        } else {
                            if ((float) (*ctx_avrQRSRR) / (float) MultiQrstFeatures[peakLoc].rr < 1.8 &&
                                (float) (*ctx_avrQRSRR) / (float) MultiQrstFeatures[peakLoc].rr > 1.0) {
                                inputData = (float) (MultiQrstFeatures[peakNumberPre].rr - (*ctx_avrQRSRR)) /
                                            (float) (*ctx_avrQRSRR);
                            } else {
                                if ((*ctx_complementFlag) == 3) {
                                    inputData =
                                            (float) (MultiQrstFeatures[peakNumberPre].rr - (*ctx_avrQRSRR)) /
                                            (float) (*ctx_avrQRSRR);
                                } else {
                                    peakTempBuffer = ((*ctx_MultiPeakNumber) - 2);
                                    if (peakTempBuffer < 0) {
                                        peakTempBuffer = peakTempBuffer + maxQueueLength;
                                    }//end

                                    inputData = (float) (MultiQrstFeatures[peakNumberPre].rr -
                                                         MultiQrstFeatures[peakLoc].rr) /
                                                (float) MultiQrstFeatures[peakTempBuffer].rr;
                                }//end
                            }
                        }//end   if(((*ctx_complementFlag)==2)&&(((rrTempInterval>0.9)&&(rrTempInterval<1.17))||rrTempInterval<0.88)||rrTempInterval>1.2)
                    }//end if(((*ctx_MultiPeakNumber)>=3)&&((*ctx_complementFlag)==0)&&((MultiQrstFeatures[peakLoc].rr)>(*ctx_avrQRSRR) ))//%%%  (*ctx_avrQRSRR) (*ctx_complementFlag)
                    //	 };//end

//  % end %% end if >=3

                    if ((inputData >= -1) && (inputData <= 3)) {
                        qrsIntervalStrong = (float) FuzzyReasoning(3, inputData);
                        inputTemp = (int) inputData;
                        qrsIntervalWeak = (float) FuzzyReasoning(4, inputData);
                    }//end

                    inputData =
                            (float) MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr /
                            (float) MultiQrstFeatures[peakNumberPre].rr;
                    //% complementData=inputData;
                    if ((inputData >= 0) && (inputData <= 7)) {
                        qrsComplementWeak = (float) FuzzyReasoning(5, inputData);
                        qrsComplementStrong = (float) FuzzyReasoning(6, inputData);
                    }//end

                    fuzzyNormalResult = (float) (qrsWidthNarrow * qrsIntervalWeak * qrsComplementWeak);
                    fuzzyPVCResult = (float) (qrsWidthWide * qrsIntervalStrong * qrsComplementStrong);
                    float rrPVCLogic = fuzzyPVCResult;

                    areaIndex = (MultiQrstFeatures[peakNumberPre].area) / (*ctx_currentArea);
                    fuzzyApcResult = (float) (qrsWidthNarrow * qrsIntervalStrong * qrsComplementStrong);

                    float peakNegativeLogic = (float) (qrsIntervalWeak * qrsComplementWeak);
                    //%            if(((fuzzyNormalResult>0.75)||(peakNegativeLogic>0.9))&&(areaIndex<1.35))   %%%&&  0.5
                    if ((fuzzyNormalResult > 0.75) && (areaIndex < 1.35) && (areaIndex > 0.65)) {
                        //%% if((fuzzyNormalResult>0.75)&&(areaIndex<1.35))
                        // float Interval=test[peakNumberPre];//%   peakLocationMod        %%%%%%%%%%%%%%%%%%  14:04 10/12 2018 ypz
                        float start = (float) Interval - (float) SAMPLE_RATE / BEAT_SAMPLE_RATE * FIDMARK;
                        if (start < 0)
                            start = start + ECG_BUFFER_LENGTH;
                        //end
                        float finish = start + (float) SAMPLE_RATE / BEAT_SAMPLE_RATE * BEATLGTH;
                        if (finish > ECG_BUFFER_LENGTH)
                            finish = finish - ECG_BUFFER_LENGTH;
                        // end
                        // %   lastBeat = ecgCopy(ECGBuffer(2),start,finish) ;   //  ECGBufferMulti
                        int channelNo = 1;//2;
                        // MultiEcgCopy(ECGBufferMulti, lastBeat, channelNo, (int) start, (int) finish);

                        // normalTemplate= templateCopy(lastBeat ,MultiQrstFeatures[peakNumberPre].onset,MultiQrstFeatures[peakNumberPre].offset);
                        templateCopy(lastBeat, normalTemplate);
                        (*ctx_templateFlag) = 1;
                        // normalP=[];
                        // % normalP=arrayCopy(lastBeat ,qrsFeatures((*ctx_peakNumber)-1).Ponset,qrsFeatures((*ctx_peakNumber)-1).Poffset);

                    }//end
// %%%%%%    20221203  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    float recentBeat[100];
                    memset(recentBeat, 0, sizeof(float) * 100);
                    // recentBeat=zeros(1,100); memset(Pheight,0,sizeof(float)*templateLength);
                    // recentBeat = templateCopy(lastBeat ,MultiQrstFeatures[peakNumberPre].onset,MultiQrstFeatures[peakNumberPre].offset);
                    templateCopy(lastBeat, recentBeat);
                    float morphologyIndex = 0;
                    int pMorphologyIndex = 0;
                    float correlation;

                    float temp[2 * BEAT_MS50];
                    memset(temp, 0, 2 * BEAT_MS50 * sizeof(float));

                    if ((*ctx_templateFlag) > 0) {
                        float mCorrelation = correlationCoefficient(normalTemplate, recentBeat, 4 * BEAT_MS50);

                        morphologyIndex = mCorrelation;
                        int pPeakPositionTemp = MultiQrstFeatures[peakNumberPre].pPeakPositionFirst;
                        if (pPeakPositionTemp - BEAT_MS50 > 0) {
                            // temp=zeros(1,2*BEAT_MS50);
                            int pIndex = pPeakPositionTemp - BEAT_MS50;
                            float temp[500];
                            int jj = 0;
                            while (pIndex < (pPeakPositionTemp + BEAT_MS50)) {
                                temp[jj] = lastBeat[pIndex];
                                pIndex = pIndex + 1;
                                jj = jj + 1;
                            }//end
                            // pMorphologyIndex=correlationCoefficient( pwaveTemplate,temp ,2*BEAT_MS50);
                            //correlation = pWaveCorrelation( pwaveTemplate,temp,(*ctx_currentPeak) );
                        }//end
                    } else {
                        //normalTemplate = templateCopy(tempBeat ,MultiQrstFeatures((*ctx_MultiPeakNumber)).onset,MultiQrstFeatures((*ctx_MultiPeakNumber)).offset);
                        templateCopy(tempBeat, normalTemplate);
                        correlation = correlationCoefficient(normalTemplate, recentBeat, 4 * BEAT_MS50);
                        morphologyIndex = correlation;
                    }
                    // end


                    // %  1.(qrsFeatures(peakNumberPre).Pwidth)<BEAT_MS140)
                    // %  2. (((qrsFeatures(peakNumberPre).pPeak)/amplifierGainConst)<0.3)&&(((qrsFeatures(peakNumberPre).pPeak)/amplifierGainConst)>0.05)
                    // %  3. P direction  (qrsFeatures(peakNumberPre).pPeak)>0
                    // %  4.   MORPHOLOGY

                    int sinusPrInterval = MultiQrstFeatures[peakNumberPre].prInterval < BEAT_MS220;
                    int sinusRrInterval = MultiQrstFeatures[peakNumberPre].rr > BEAT_MS600
                                          && MultiQrstFeatures[peakNumberPre].rr < BEAT_MS1200;
                    float pArea = MultiQrstFeatures[peakNumberPre].Psum / (*ctx_currentPsum);
                    float pPeakIndex = MultiQrstFeatures[peakNumberPre].pPeak / (*ctx_currentPeak);
                    int sinusIndex = (pPeakIndex > 0.65 && pPeakIndex < 1.6)
                                     || (pArea > 0.63 && pArea < 1.6);

// %    qrsAreaRation=qrsFeatures(peakNumberPre).area/(*ctx_currentArea);
                    int sinusPwidth = MultiQrstFeatures[peakNumberPre].Pwidth <= BEAT_MS140;
                    int sinusAmplifierRange = MultiQrstFeatures[peakNumberPre].pPeak / amplifierGainConst < 0.3
                                              && MultiQrstFeatures[peakNumberPre].pPeak / amplifierGainConst > 0.05;
                    int sinusParameter = sinusIndex && sinusPwidth && sinusAmplifierRange;
                    int pPeakPositionTemp;
                    int pIndex;
                    if (sinusPrInterval && sinusRrInterval && sinusParameter) {
                        // %%%%%%%% p wave template  %%%%%%%%%%%%%%%%
                        pPeakPositionTemp = MultiQrstFeatures[peakNumberPre].pPeakPositionFirst;
                        if (MultiQrstFeatures[peakNumberPre].pType == 1) {
                            pIndex = pPeakPositionTemp - BEAT_MS50;
                            if (pIndex > 0) {
                                int jj = 0;
                                while (pIndex < (pPeakPositionTemp + BEAT_MS50)) {
                                    pwaveTemplate[jj] = lastBeat[pIndex];
                                    pIndex = pIndex + 1;
                                    jj = jj + 1;
                                }//end
                            }//end
                        }//end
                        // %%%%%%%% p wave template end  %%%%%%%%%%%%%%%%
                    }//end

                    //     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ectopicIndex  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    pArea = (MultiQrstFeatures[peakNumberPre].Psum / (*ctx_currentPsum));
                    pPeakIndex = (MultiQrstFeatures[peakNumberPre].pPeak) / (*ctx_currentPeak);
                    sinusIndex = ((pPeakIndex > 0.85) && (pPeakIndex < 1.15)) || ((pArea > 0.85) && (pArea < 1.15));
                    float qrsAreaRation = MultiQrstFeatures[peakNumberPre].area / (*ctx_currentArea);

                    //%  qrsParameter=((abs(((qrsFeatures(peakNumberPre).pPeak))/(*ctx_currentPeak))<6.8)&&((qrsFeatures(peakNumberPre).qrsWidth)<130)&&(qrsFeatures(peakNumberPre).area)<1.5*ctx_maxQRSAerea) ;
                    int qrsParameter = MultiQrstFeatures[peakNumberPre].qrsWidth < 130
                                       &&
                                       MultiQrstFeatures[peakNumberPre].area <
                                       1.5 * (*ctx_maxQRSAerea);// %%%  20191216  ypz
                    int ectopicIndex = MultiQrstFeatures[peakNumberPre].pPeak / (*ctx_currentPeak) > 1.2
                                       && MultiQrstFeatures[peakNumberPre].Psum / (*ctx_currentPsum) > 1.28
                                       || MultiQrstFeatures[peakNumberPre].pPeak / (*ctx_currentPeak) > 1.7;

                    int tWaveDirection =
                            MultiQrstFeatures[peakNumberPre].tHeight * MultiQrstFeatures[peakNumberPre].qrsMainWave < 0;
                    int ventricularIndex =
                            fabsf(MultiQrstFeatures[peakNumberPre].area) > 1.8 * (*ctx_currentArea) &&
                            sinusIndex == 0 &&
                            MultiQrstFeatures[peakNumberPre].rr < BEAT_MS800 &&
                            (MultiQrstFeatures[peakNumberPre].tAera / (*ctx_currentTArea) > 1.65 ||
                             fabsf(MultiQrstFeatures[peakNumberPre].speak) > 1.55 * (*ctx_currentAverSpeak));

                    //%  ventricularIndex=(((abs(qrsFeatures(peakNumberPre).area))>1.8*ctx_currentArea)&&(sinusIndex==0)&&((qrsFeatures(peakNumberPre).rr)<BEAT_MS800)&&((((qrsFeatures(peakNumberPre).tAera)/(*ctx_currentTArea))>1.65)||(abs(qrsFeatures(peakNumberPre).speak))>1.55*ctx_currentAverSpeak)&&tWaveDirection);
                    //%
                    float h = MultiQrstFeatures[peakNumberPre].area / (*ctx_currentArea);
                    float g = MultiQrstFeatures[peakNumberPre].tAera / (*ctx_currentTArea);
                    float k = fabsf(MultiQrstFeatures[peakNumberPre].speak) / (*ctx_currentAverSpeak);

                    //%  if((qrsFeatures(peakNumberPre).pType==2)&&(qrsFeatures(peakNumberPre).prInterval<BEAT_MS120)&&(morphologyIndex>=0.92)&&(qrsIntervalStrong>0.6)&&((*ctx_currentPeak)/abs(qrsFeatures(peakNumberPre).pPeak))>1.5)
                    /*  if((MultiQrstFeatures[peakNumberPre].pType==2)&&(MultiQrstFeatures[peakNumberPre].prInterval<BEAT_MS120)&&(morphologyIndex>=0.92)&&((*ctx_currentPeak)/abs(MultiQrstFeatures[peakNumberPre].pPeak))>1.5){
       jPremature=jPremature+1;
       japcFlag=1;
       [ jpcBigeminyCount,jpcTrigeminyCount,jpcSupraventricularCount] = junctionalPrematureCount( japcFlag  ,0);

   }//end

       if((MultiQrstFeatures(peakNumberPre).pType==2)&&(MultiQrstFeatures(peakNumberPre).prInterval<BEAT_MS120)&&(morphologyIndex>=0.92)&&((*ctx_currentPeak)/abs(MultiQrstFeatures(peakNumberPre).pPeak))>1.5)
          jPremature=jPremature+1;
          japcFlag=1;
          [ jpcBigeminyCount,jpcTrigeminyCount,jpcSupraventricularCount] = junctionalPrematureCount( japcFlag  ,0);

       end
*/
//%  hugeVentricular=((qrsFeatures(peakNumberPre).area)>3.3*ctx_currentArea)&&((qrsFeatures(peakNumberPre).pType==5)||(qrsFeatures(peakNumberPre).pType==0))&&(inputTemp<-0.05)&&(inputTemp>-0.2);
//%%hugeVentricular=((qrsFeatures(peakNumberPre).area)>3.3*ctx_currentArea)&&((qrsFeatures(peakNumberPre).pType==5)||((qrsFeatures(peakNumberPre).qrsWidth)>BEAT_MS160)&&((((qrsFeatures(peakNumberPre).tAera)/(*ctx_currentTArea))>1.65)&&(abs(qrsFeatures(peakNumberPre).speak))>1.55*ctx_currentAverSpeak)||(qrsFeatures(peakNumberPre).pType==0))&&(inputTemp<0.15)&&(inputTemp>-0.2);
                    int hugeVentricular = MultiQrstFeatures[peakNumberPre].area > 3.0 * (*ctx_currentArea)
                                          && (MultiQrstFeatures[peakNumberPre].pType == 5
                                              || (MultiQrstFeatures[peakNumberPre].qrsWidth > BEAT_MS160
                                                  &&
                                                  (MultiQrstFeatures[peakNumberPre].tAera / (*ctx_currentTArea) > 1.65
                                                   || fabsf(MultiQrstFeatures[peakNumberPre].speak) >
                                                      1.55 * (*ctx_currentAverSpeak)))
                                              || MultiQrstFeatures[peakNumberPre].pType == 0)
                                          && inputTemp < 0.15
                                          && inputTemp > -0.2;
                    if (hugeVentricular == 1)
                        ventricularIndex = 1;
//                        ventricularIndex = 1;
//end
                    int index = 0;
                    float normalIndex = 0;
                    while (index < apcProcessLength) {
                        apcProcessing[index] = apcProcessing[index + 1];
                        index = index + 1;
                    }//end
                    normalIndex = normalIndex / apcProcessLength;
                    float RRav = (float) MultiQrstFeatures[peakNumberPre].rr / SAMPLE_RATE;
                    float HR = 60 / RRav;
                    if (((*ctx_complementFlag) == 2) && (peakNumberPre > 34)) {
                        apcProcessing[apcProcessLength - 1] = 0;
                        apcProcessing[apcProcessLength] = 0;
                    } else {
                        if (((MultiQrstFeatures[peakNumberPre].rr) <= BEAT_MS1600))
                            apcProcessing[apcProcessLength] = (float) MultiQrstFeatures[peakNumberPre].rr;
                        else
                            apcProcessing[apcProcessLength] = 0;
                        //end
                    }//end

//% [shannonEntropy] =  shannonDistribution(apcProcessing);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ectopic
                    int apcFlag = 0;
                    float ll = (float) MultiQrstFeatures[peakNumberPre].qrsWidth;

                    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    int Fpa = 50;
                    //int Fs=500;
                    float X[SAMPLE_RATE];
                    //float X  = LowPassFilter(lastBeat,Fs,Fpa);
                    LowPassFilter( lastBeat, X);
                    int start = 60;
                    int finish = 155;//%165;
                    memset(newArray, 0, BEAT_MS200 * sizeof(float));
                    ecgCopy(X, newArray, start, finish);
                    int len = finish - start;
                    float mnoise = noiseLevel(newArray, len);//letHeight = {float} -24.1806641n);
                    MultiQrstFeatures[peakNumberPre].noiseCount = (int) mnoise;
                    if ((mnoise >= 20) && (MultiQrstFeatures[peakNumberPre].qrsHeight > 0)) {
                        return 0;
                        // %continue;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        //%
                    }//end
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    // int apcFlag=0;
                    int japcFlag = 0;//  !!!!!!!!!!!!!!!!!!!!!!!!

                    int temppp = peakNumberPre;
                    if (temppp == 305) {
                        temppp = 0;
                    }
                    if (((qrsIntervalStrong > 0.749) || ((ventricularIndex == 1) && (qrsParameter == 1))) &&
                        (MultiQrstFeatures[peakNumberPre].qrsWidth > 0)) {
                        shannonEntropy = shannonDistribution( apcProcessing, apcProcessLength);
                        if (((*ctx_currentRRSMD) > 0.1) && (shannonEntropy < (*ctx_currentRRSMD))) {
                            shannonEntropy = RMSSD;
                        }
                        float pvcTemplateLength = 130;
                        int templateCenter = 140;
                        int m;
                        if (((*ctx_pvcLength) > pvcQueLength) && (((*ctx_complementFlag) != 2)) &&
                            (MultiQrstFeatures[peakNumberPre].area > (float) (*ctx_currentWidthMax)) && (ll >= 75) &&
                            (((ventricularIndex == 1)) ||
                             (MultiQrstFeatures[peakNumberPre].qrsWidth > (*ctx_currentQRS)))) {
                            m = 0;
                            while ((float) m <= pvcTemplateLength) {
                                newArray[m] = lastBeat[m + templateCenter];
                                m = m + 1;
                            }//end
                            float maxCoef = templateSeek( &(*ctx_pvcLibInfo), newArray);
                            if ((maxCoef > 0.92))//%0.89))
                            {
                                (*ctx_pvcNumber) = (*ctx_pvcNumber) + 1;
                                (*ctx_complementFlag) = 1;
                                pvcFlag = 1;
                            }//end
                            continue;
                        }//end

                        if (((*ctx_complementFlag) == 3) && (qrsWidthWide == 1) && (qrsIntervalStrong == 1) &&
                            (qrsComplementStrong < 0.6))//%% for short complement period
                            fuzzyPVCResult = 0.93f;
                        //end
                        float ectopic = (MultiQrstFeatures[peakNumberPre].pPeak) / (*ctx_tHeightAverage);

                        int pvcAfFlag = 1;
                        if ((shannonEntropy > 0.1) && ((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS120))
                            pvcAfFlag = 0;
                        //end



                        //%%% 11/12   2022      ypz

                        if ((MultiQrstFeatures[peakLoc].noiseFlag == 0) &&
                            (((MultiQrstFeatures[peakNumberPre].pPeak) / (MultiQrstFeatures[peakLoc].tHeight)) >
                             0.95) &&
                            (((MultiQrstFeatures[peakNumberPre].pPeak) / (MultiQrstFeatures[peakLoc].tHeight)) < 1.05))
                            //%       ptemp=(((qrsFeatures((*ctx_peakNumber)).pPeak))/amplifierGain>0.4);%gainRate=0.4
                            //  %  if(ptemp&&(qrsFeatures((*ctx_peakNumber)-2).tFirstPeak>1)&&(qrsFeatures((*ctx_peakNumber)-1).pPeakPositionFirst>1))
                        {
                            // Fusion beat between two beats(T Wave and P Wave )
                            int len;
                            if ((MultiQrstFeatures[peakLoc].tFirstPeak > 1) &&
                                (MultiQrstFeatures[peakNumberPre].pPeakPositionFirst > 1)) {
                                len = 25;
                                float correlation = arrayCreate(secondBeat, lastBeat,
                                                                MultiQrstFeatures[peakLoc].tFirstPeak,
                                                                MultiQrstFeatures[peakNumberPre].pPeakPositionFirst,
                                                                len);
                                // % areaIndex=(qrsFeatures((*ctx_peakNumber)-1).area)/(*ctx_currentArea);
                                if (correlation > 0.96) {
                                    // %   areaIndex=(qrsFeatures(peakNumberPre).area)/(*ctx_currentArea);
                                    //%    ectopic=(qrsFeatures((*ctx_peakNumber)-1).pPeak)/(*ctx_tHeightAverage);%(qrsFeatures((*ctx_peakNumber)-1).area))>1.40*ctx_currentArea)
                                    //%   if((ectopic>1.2)&&((qrsFeatures((*ctx_peakNumber)-1).qrsWidth)<55)&&( morphologyIndex>0.95))
                                    if ((((ectopic > 1.2) && (shannonEntropy < 0.1) &&
                                          ((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS110)) ||
                                         (MultiQrstFeatures[peakNumberPre].pType == 5) || (ectopicIndex == 1) ||
                                         (MultiQrstFeatures[peakNumberPre].pType == 6) ||
                                         ((MultiQrstFeatures[peakNumberPre].prInterval) >= 100)) &&
                                        (morphologyIndex > 0.92) && (areaIndex < 1.31)) {
                                        (*ctx_apc) = (*ctx_apc) + 1;
                                        apcFlag = 1;
                                        (*ctx_complementFlag) = 1;
                                    } else {
                                        //                   if(areaIndex>1.4)&&(ectopic>1.2)&&( morphologyIndex<0.95)
                                        //                    % aberrantConduction=1;
//%                      (*ctx_apc)=(*ctx_apc)+1;
//%                      apcFlag=1;
//%                      (*ctx_complementFlag)=1;
//                    end

                                        if ((ectopic > 1.2) && (shannonEntropy < 0.1) && (areaIndex < 1.5) &&
                                            ((MultiQrstFeatures[peakNumberPre].qrsWidth) <= BEAT_MS110) &&
                                            (morphologyIndex > 0.92))//%%% &&( morphologyIndex>0.95)
                                        {
                                            (*ctx_apc) = (*ctx_apc) + 1;
                                            apcFlag = 1;
                                            (*ctx_complementFlag) = 1;
                                        } else {

                                            float c, d;
                                            if (((*ctx_pvcSum) >= 4)) {
                                                c = (MultiQrstFeatures[peakNumberPre].area) / (*ctx_pvcStandard);
                                                d = (float) abs(MultiQrstFeatures[peakNumberPre].qrsWidth) /
                                                    (*ctx_pvcWidthAverage);
                                                // %  k=(c*d<1.9)&&(c*d>0.6);

                                            }// end

                                            if ((((areaIndex > 1.5) || ((fuzzyPVCResult > 0.70))) && (pvcAfFlag == 1) &&
                                                 (morphologyIndex < 0.96)) &&
                                                (((*ctx_pvcSum) >= 4) && ((c * d < 1.9) || ((*ctx_pvcSum) < 4))))
                                                // %% %  if(((areaIndex>1.5)||((fuzzyPVCResult>0.70)))&&( morphologyIndex<0.92))
                                            {
                                                (*ctx_pvc) = (*ctx_pvc) + 1;
                                                (*ctx_complementFlag) = 1;
                                                pvcFlag = 1;

                                                if ((*ctx_pvcSum) < 4) {
                                                    (*ctx_pvcSum) = (*ctx_pvcSum) + 1;
                                                }//end

                                                index = 0;
                                                while (index < 4) {
                                                    pvcArray[index].qrsArea = pvcArray[index + 1].qrsArea;
                                                    pvcArray[index].qrsWidth = pvcArray[index + 1].qrsWidth;
                                                    index = index + 1;
                                                }//end

                                                pvcArray[3].qrsArea = (MultiQrstFeatures[peakNumberPre].area);
                                                pvcArray[3].qrsWidth = (MultiQrstFeatures[peakNumberPre].qrsWidth);

                                                if ((*ctx_pvcSum) >= 3) {
                                                    (*ctx_pvcWidthAverage) =
                                                            (float) (pvcArray[0].qrsWidth + pvcArray[1].qrsWidth +
                                                                     pvcArray[2].qrsWidth + pvcArray[3].qrsWidth) / 4;
                                                    (*ctx_pvcStandard) = (pvcArray[0].qrsArea + pvcArray[1].qrsArea +
                                                                          pvcArray[2].qrsArea + pvcArray[3].qrsArea) /
                                                                         4;
                                                    (*ctx_pvcStandardDerivation) = sqrtf(
                                                            (
                                                                    (pvcArray[0].qrsArea - (*ctx_pvcStandard)) *
                                                                    (pvcArray[0].qrsArea - (*ctx_pvcStandard))
                                                                    +
                                                                    (pvcArray[1].qrsArea - (*ctx_pvcStandard)) *
                                                                    (pvcArray[1].qrsArea - (*ctx_pvcStandard))
                                                                    +
                                                                    (pvcArray[2].qrsArea - (*ctx_pvcStandard)) *
                                                                    (pvcArray[2].qrsArea - (*ctx_pvcStandard))
                                                                    +
                                                                    (pvcArray[3].qrsArea - (*ctx_pvcStandard)) *
                                                                    (pvcArray[3].qrsArea - (*ctx_pvcStandard))
                                                            ) / 4);
                                                    (*ctx_pvcStandardWidth) = sqrtf(
                                                            (
                                                                    ((float) pvcArray[0].qrsWidth -
                                                                     (*ctx_pvcWidthAverage)) *
                                                                    ((float) pvcArray[0].qrsWidth -
                                                                     (*ctx_pvcWidthAverage))
                                                                    +
                                                                    ((float) pvcArray[1].qrsWidth -
                                                                     (*ctx_pvcWidthAverage)) *
                                                                    ((float) pvcArray[1].qrsWidth -
                                                                     (*ctx_pvcWidthAverage))
                                                                    +
                                                                    ((float) pvcArray[2].qrsWidth -
                                                                     (*ctx_pvcWidthAverage)) *
                                                                    ((float) pvcArray[2].qrsWidth -
                                                                     (*ctx_pvcWidthAverage))
                                                                    +
                                                                    ((float) pvcArray[3].qrsWidth -
                                                                     (*ctx_pvcWidthAverage)) *
                                                                    ((float) pvcArray[3].qrsWidth -
                                                                     (*ctx_pvcWidthAverage))
                                                            ) / 4);
                                                    if ((*ctx_pvcStandardWidth) > (*ctx_pvcWidthMaximum))
                                                        (*ctx_pvcWidthMaximum) = (*ctx_pvcStandardWidth);
                                                    //end
                                                    if ((*ctx_pvcStandardDerivation) > (*ctx_pvcAreaMaximum))
                                                        (*ctx_pvcAreaMaximum) = (*ctx_pvcStandardDerivation);
                                                    //end
                                                }//end

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                int m = 0;
                                                //int (*ctx_n)=0;
                                                int k;
                                                if ((*ctx_n) <= 1) {
                                                    while (m < 100) {
                                                        (*ctx_pvcLibInfo).informArray[1][m] = lastBeat[m + 170];
                                                        m = m + 1;
                                                    }//end
                                                    (*ctx_pvcLibInfo).queueLength = 1;
                                                }//end
                                                (*ctx_n) = (*ctx_n) + 1;

                                                if ((*ctx_n) > 2) {
                                                    k = 0;
                                                    while (k < 100) {
                                                        newArray[k] = lastBeat[k + 170];
                                                        k = k + 1;
                                                    }//end
                                                    templateUpdate( &(*ctx_pvcLibInfo), newArray);
                                                }//end

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            }//    end  if(((areaIndex>1.5)||((fuzzyPVCResult>0.70)))&&(pvcAfFlag==1)&&( morphologyIndex<0.96))&&(((*ctx_pvcSum)>=4)&&((c*d<1.9)||((*ctx_pvcSum)<4)))
                                        }// end  4462          else
                                    }// line 4447 end   if(((ectopic>1.2)&&(shannonEntropy<0.1)&&((MultiQrstFeatures[peakNumberPre].

                                }//     end   if(correlation>0.96)
                                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if ((apcFlag == 0) && (pvcFlag == 0) &&
                                    ((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS150)) {
                                    float qrsAreaRation = MultiQrstFeatures[peakNumberPre].area / (*ctx_currentArea);
                                    int apcLogicOne = ((((qrsAreaRation < 1.1)) &&
                                                        ((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS110)) ||
                                                       ((ectopicIndex == 1) && (morphologyIndex > 0.92))) &&
                                                      (qrsIntervalStrong == 1);
                                    int pWaveFeature = ((MultiQrstFeatures[peakNumberPre].Psum / (*ctx_currentPsum)) &&
                                                        ((MultiQrstFeatures[peakNumberPre].pPeak) > 0));
                                    int apcLogicTwo = (((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS120) &&
                                                       ((MultiQrstFeatures[peakNumberPre].pType == 5) ||
                                                        (pWaveFeature > 2)) && (morphologyIndex > 0.87));
                                    if (((apcLogicOne || apcLogicTwo) && (ventricularIndex == 0)) &&
                                        (shannonEntropy < 0.1)) {
                                        (*ctx_apc) = (*ctx_apc) + 1;
                                        apcFlag = 1;
                                        (*ctx_complementFlag) = 1;
                                        // % if(qrsFeatures(peakNumberPre).area>(*ctx_maxQRSAerea))
                                        //  %   (*ctx_maxQRSAerea)=qrsFeatures(peakNumberPre).area;
                                        //  %end
                                    }//end if((apcLogicOne||apcLogicTwo)&&(ventricularIndex==0))&&(shannonEntropy<0.1)
                                }//end if((apcFlag==0)&&(pvcFlag==0)&&((MultiQrstFeatures(peakNumberPre).qrsWidth)<BEAT_MS150))
                                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            }// 4432 end  if((MultiQrstFeatures[peakLoc].tFirstPeak>1)
                            float c, d;
                            if (((*ctx_pvcSum) >= 4)) {
                                c = (MultiQrstFeatures[peakNumberPre].area) / (*ctx_pvcStandard);
                                d = (float) abs((MultiQrstFeatures[peakNumberPre].qrsWidth)) / (*ctx_pvcWidthAverage);
                            }//end

                            if ((((ventricularIndex == 1) || ((morphologyIndex < 0.93) && (fuzzyPVCResult > 0.5) &&
                                                              ((qrsIntervalStrong) == 1.0))) && (pvcFlag == 0) &&
                                 (pvcAfFlag == 1) && (apcFlag == 0)) && ((((*ctx_pvcSum) >= 4) /*&& (c * d < 1.9)*/) ||
                                                                         ((morphologyIndex < 0.93) &&
                                                                          (fuzzyPVCResult == 1)) ||
                                                                         ((*ctx_pvcSum) < 4)))
                                // % if(((ventricularIndex==1)||(morphologyIndex<0.9)&&(fuzzyPVCResult>0.5)&&((qrsIntervalStrong)==1.0))&&(pvcFlag==0)&&(apcFlag==0))
                            {
                                (*ctx_pvc) = (*ctx_pvc) + 1;
                                (*ctx_complementFlag) = 1;
                                pvcFlag = 1;

                                //    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                int m = 0, k;
                                if ((*ctx_n) < 1) {
                                    while (m < 100) {
                                        (*ctx_pvcLibInfo).informArray[1][m] = lastBeat[m + 170];
                                        m = m + 1;
                                    }// end
                                    (*ctx_pvcLibInfo).queueLength = 1;
                                }//end
                                (*ctx_n) = (*ctx_n) + 1;

                                if ((*ctx_n) > 2) {
                                    k = 1;
                                    while (k < 100) {
                                        newArray[k] = lastBeat[k + 170];
                                        k = k + 1;
                                    }//end
                                    // [ (*ctx_pvcLibInfo) ] = templateUpdate( (*ctx_pvcLibInfo),newArray );
                                    templateUpdate(&(*ctx_pvcLibInfo), newArray);
                                }//end

                                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                if ((*ctx_pvcSum) < 4)
                                    (*ctx_pvcSum) = (*ctx_pvcSum) + 1;
                                // end

                                index = 0;
                                while (index < 4) {
                                    pvcArray[index].qrsArea = pvcArray[index + 1].qrsArea;
                                    pvcArray[index].qrsWidth = pvcArray[index + 1].qrsWidth;
                                    index = index + 1;
                                }//end

                                pvcArray[3].qrsArea = (MultiQrstFeatures[peakNumberPre].area);
                                pvcArray[3].qrsWidth = (MultiQrstFeatures[peakNumberPre].qrsWidth);

                                if ((*ctx_pvcSum) >= 3) {
                                    (*ctx_pvcWidthAverage) =
                                            (float) (pvcArray[0].qrsWidth + pvcArray[1].qrsWidth +
                                                     pvcArray[2].qrsWidth +
                                                     pvcArray[3].qrsWidth) / 4;
                                    (*ctx_pvcStandard) =
                                            (pvcArray[0].qrsArea + pvcArray[1].qrsArea + pvcArray[2].qrsArea +
                                             pvcArray[3].qrsArea) / 4;
                                    (*ctx_pvcStandardDerivation) = sqrtf(
                                            ((pvcArray[0].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[0].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[1].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[1].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[2].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[2].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[3].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[3].qrsArea - (*ctx_pvcStandard))) / 4);
                                    (*ctx_pvcStandardWidth) = sqrtf(
                                            (((float) pvcArray[0].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[0].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[1].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[1].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[2].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[2].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[3].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[3].qrsWidth - (*ctx_pvcWidthAverage))) / 4);
                                    if ((*ctx_pvcStandardWidth) > (*ctx_pvcWidthMaximum))
                                        (*ctx_pvcWidthMaximum) = (*ctx_pvcStandardWidth);
                                    //end
                                    if ((*ctx_pvcStandardDerivation) > (*ctx_pvcAreaMaximum))
                                        (*ctx_pvcAreaMaximum) = (*ctx_pvcStandardDerivation);
                                    //end
                                }//end


                            }//end    if(((ventricularIndex==1)||(morphologyIndex<0.932)&&(fuzzyPVCResult>0.5)&&((qrsIntervalStrong)==1.0))&
                        }//end      %%%   if((MultiQrstFeatures(peakLoc).noiseFlag==0)



                        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& PVC &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        if ((pvcAfFlag == 1) && (apcFlag == 0) && (pvcFlag == 0) && ((morphologyIndex < 0.93) ||
                                                                                     ((fuzzyPVCResult >= 0.98) &&
                                                                                      (morphologyIndex <
                                                                                       0.94))))//%%&&((qrsFeatures(peakNumberPre).area)<1.5*ctx_maxQRSAerea))
                            //% areaIndex=(qrsFeatures(peakNumberPre).area)/(*ctx_currentArea);
                            //  float a,b,c,d;
                        {
                            float c = 0;
                            float d;
                            int k;
                            // a=abs((MultiQrstFeatures[peakNumberPre].area)-(*ctx_pvcStandardDerivation))/(*ctx_pvcAreaMaximum);
                            // b=(abs((MultiQrstFeatures[peakNumberPre].qrsWidth)-(*ctx_pvcWidthAverage))/(*ctx_pvcWidthMaximum));
                            if (((*ctx_pvcSum) >= 4)) {
                                c = (MultiQrstFeatures[peakNumberPre].area) / (*ctx_pvcStandard);
                                d = (float) abs(MultiQrstFeatures[peakNumberPre].qrsWidth) / (*ctx_pvcWidthAverage);
                                k = c * d < 1.9 && c * d > 0.6;


                            }//end

                            // %  pvcTag=(((qrsWidthWide*qrsIntervalStrong)>=0.73)&&(abs(qrsFeatures(peakNumberPre).pPeak)<5.0)&&(((qrsFeatures(peakNumberPre).tAera)/(*ctx_currentTArea))>2.5)&&(((qrsFeatures(peakNumberPre).speak)/(*ctx_currentAverSpeak))>2.5||(qrsFeatures(peakNumberPre).qrsWidth)>=BEAT_MS140));
                            int pvcNoiseFlag = (qrsIntervalStrong < 0.1) &&
                                               ((fabsf((MultiQrstFeatures[peakNumberPre].area) -
                                                       (*ctx_pvcStandardDerivation)) / (*ctx_pvcAreaMaximum)) > 1.5) &&
                                               ((fabsf((float) MultiQrstFeatures[peakNumberPre].qrsWidth -
                                                       (*ctx_pvcWidthAverage)) /
                                                 (*ctx_pvcWidthMaximum)) > 1.5);
//%      if(((morphologyIndex<0.94)||(fuzzyPVCResult>=0.8)&&(areaIndex>1.15)))&&(((*ctx_pvcSum)>=4)&&(c*d<1.9)||((*ctx_pvcSum)<4))
                            int pvcTag = (((qrsWidthWide * qrsIntervalStrong) >= 0.725) &&
                                          (fabsf(MultiQrstFeatures[peakNumberPre].pPeak) < 5.0) &&
                                          (((MultiQrstFeatures[peakNumberPre].tAera) / (*ctx_currentTArea)) > 2.5) &&
                                          (((MultiQrstFeatures[peakNumberPre].speak) / (*ctx_currentAverSpeak)) > 2.5 ||
                                           ((morphologyIndex < 0) && (tWaveDirection == 1)) ||
                                           (MultiQrstFeatures[peakNumberPre].qrsWidth) >= BEAT_MS140));
                            pvcFlag = (((((qrsWidthWide * qrsIntervalStrong) > 0.85) &&
                                         (MultiQrstFeatures[peakNumberPre].pType == 0) && (qrsAreaRation > 1.10)) ||
                                        (MultiQrstFeatures[peakNumberPre].area > 2.0 * (*ctx_currentArea)) ||
                                        (MultiQrstFeatures[peakNumberPre].pType == 5)));//%
                            if ((fuzzyPVCResult > 0.85) || pvcTag || pvcFlag || (ventricularIndex == 1) ||
                                ((fuzzyPVCResult > 0.5) && (((((MultiQrstFeatures[peakNumberPre].qrsMainWave) *
                                                               (MultiQrstFeatures[peakNumberPre].tHeight)) < 0) &&
                                                             ((MultiQrstFeatures[peakNumberPre].tAera) /
                                                              (*ctx_currentTArea)) >
                                                             1.8) ||
                                                            ((MultiQrstFeatures[peakNumberPre].tAera) /
                                                             (*ctx_currentTArea)) >
                                                            4 || (fabsf(MultiQrstFeatures[peakNumberPre].area) >
                                                                  1.40 * (*ctx_currentArea)) ||
                                                            (((MultiQrstFeatures[peakNumberPre].qrsMainWave) > 0) &&
                                                             (fabsf(MultiQrstFeatures[peakNumberPre].speak) >
                                                              1.50 * (*ctx_currentAverSpeak))))))//%pvcThrehold)
                            {
                                if (((*ctx_sIndex) < 2) && ((*ctx_currentAverSpeak) == maxInteger)) {
                                    (*ctx_currentAverSpeak) =
                                            fabsf(MultiQrstFeatures[peakNumberPre].speak) * 4 / 7;//%*0.8/1.4;
                                    (*ctx_currentArea) = fabsf(MultiQrstFeatures[peakNumberPre].area) / 2;
                                }//end


                                if ((*ctx_maxQRSAerea) == maxInteger)
                                    (*ctx_maxQRSAerea) = MultiQrstFeatures[peakNumberPre].area;
                                else if (MultiQrstFeatures[peakNumberPre].area > (*ctx_maxQRSAerea))
                                    (*ctx_maxQRSAerea) = MultiQrstFeatures[peakNumberPre].area;
                                //end
                                //end

                                if ((*ctx_pvcSum) < 4)
                                    (*ctx_pvcSum) = (*ctx_pvcSum) + 1;
                                // end

                                index = 0;
                                while (index < 4) {
                                    pvcArray[index].qrsArea = pvcArray[index + 1].qrsArea;
                                    pvcArray[index].qrsWidth = pvcArray[index + 1].qrsWidth;
                                    index = index + 1;
                                }//end

                                pvcArray[3].qrsArea = (MultiQrstFeatures[peakNumberPre].area);
                                pvcArray[3].qrsWidth = (MultiQrstFeatures[peakNumberPre].qrsWidth);

                                if ((*ctx_pvcSum) >= 3) {
                                    (*ctx_pvcWidthAverage) =
                                            (float) (pvcArray[0].qrsWidth + pvcArray[1].qrsWidth +
                                                     pvcArray[2].qrsWidth +
                                                     pvcArray[3].qrsWidth) / 4;
                                    (*ctx_pvcStandard) =
                                            (pvcArray[0].qrsArea + pvcArray[1].qrsArea + pvcArray[2].qrsArea +
                                             pvcArray[3].qrsArea) / 4;
                                    (*ctx_pvcStandardDerivation) = sqrtf(
                                            ((pvcArray[0].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[0].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[1].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[1].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[2].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[2].qrsArea - (*ctx_pvcStandard)) +
                                             (pvcArray[3].qrsArea - (*ctx_pvcStandard)) *
                                             (pvcArray[3].qrsArea - (*ctx_pvcStandard))) / 4);
                                    (*ctx_pvcStandardWidth) = sqrtf(
                                            (((float) pvcArray[0].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[0].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[1].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[1].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[2].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[2].qrsWidth - (*ctx_pvcWidthAverage)) +
                                             ((float) pvcArray[3].qrsWidth - (*ctx_pvcWidthAverage)) *
                                             ((float) pvcArray[3].qrsWidth - (*ctx_pvcWidthAverage))) / 4);
                                    if ((*ctx_pvcStandardWidth) > (*ctx_pvcWidthMaximum))
                                        (*ctx_pvcWidthMaximum) = (*ctx_pvcStandardWidth);
                                    //end
                                    if ((*ctx_pvcStandardDerivation) > (*ctx_pvcAreaMaximum))
                                        (*ctx_pvcAreaMaximum) = (*ctx_pvcStandardDerivation);
                                    //end
                                }//end


                                int m = 0;
                                pvcFlag = 1;
                                //%    formFalg=1;
                                fuzzyPVCResult = 0;
                                (*ctx_pvcLibInfo).Sdeep[(*ctx_n)] = fabsf(MultiQrstFeatures[peakNumberPre].speak);
                                (*ctx_pvc) = (*ctx_pvc) + 1;
                                (*ctx_complementFlag) = 1;
                                if ((*ctx_n) < 1) {
                                    while (m <= 100) {
                                        (*ctx_pvcLibInfo).informArray[1][m] = lastBeat[m + 170];
                                        m = m + 1;
                                    }//end
                                    (*ctx_pvcLibInfo).queueLength = 1;
                                }//end
                                (*ctx_n) = (*ctx_n) + 1;

                                if ((*ctx_n) > 2) {
                                    int k = 0;
                                    while (k < 100) {
                                        newArray[k] = lastBeat[k + 170];
                                        k = k + 1;
                                    }//end
                                    // [ (*ctx_pvcLibInfo) ] = templateUpdate( (*ctx_pvcLibInfo),newArray );
                                    templateUpdate(&(*ctx_pvcLibInfo), newArray);
                                }// end
                            }//end

                        }//  end   %%  end  if((pvcAfFlag==1)&&(apcFlag==0)&&(pvcFlag==0)&&((morphologyIndex<0.93)

                        //%%   if((pvcFlag==0)&



                        //   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APC  %%%%%%%%%%%%%%%%%%%%%%%%%%
                        //%  (*ctx_testCount)=(*ctx_testCount)+1;
                        if ((apcFlag == 0) && (pvcFlag == 0) && (shannonEntropy < 0.1) &&
                            ((areaIndex < 1.4) || (MultiQrstFeatures[peakNumberPre].qrsWidth <= 43)) &&
                            (morphologyIndex >= 0.92) && (ventricularIndex == 0) && (japcFlag == 0)) {
                            //%       if(( sinusParameter==0)&&(pvcFlag==0)&&(shannonEntropy<0.1)&&(areaIndex<1.4)&&(morphologyIndex>=0.92)&&(ventricularIndex==0)&&(japcFlag==0))%%&&(qrsIntervalStrong==1)) %%%%%%%%%%%%  APC %%%%%%%%%%%%
                            //%        if((apcFlag==0)&&(pvcFlag==0)&&(shannonEntropy<0.1)&&((areaIndex<1.4)||(areaIndex<1.52)&&(morphologyIndex>=0.95))&&(morphologyIndex>=0.92)&&(ventricularIndex==0)&&(japcFlag==0))
                            int kk = ((*ctx_currentPeak) / MultiQrstFeatures[peakNumberPre].pPeak) > 2;
                            int m = ((MultiQrstFeatures[peakNumberPre].pPeak) / (*ctx_currentPeak) > 1.24) ||
                                    (((*ctx_currentPeak) / MultiQrstFeatures[peakNumberPre].pPeak) > 2);
                            if ((*ctx_MultiPeakNumber) <= 9) {
                                m = ((MultiQrstFeatures[peakNumberPre].pPeak) / MultiQrstFeatures[peakLoc].pPeak > 1.3);

                                peakTempBuffer = ((*ctx_MultiPeakNumber) -
                                                  2);                          // %%%%%%%%%%%%%%%%%%%
                                if (peakTempBuffer < 1) {
                                    peakTempBuffer = peakTempBuffer + maxQueueLength;
                                }//end

                                (*ctx_currentArea) = MultiQrstFeatures[peakTempBuffer].area;
                            }//end
                            int interpolatedIndex = ((MultiQrstFeatures[peakLoc].noiseFlag == 0) &&
                                                     ((float) (MultiQrstFeatures[peakNumberPre].rr +
                                                               MultiQrstFeatures[peakLoc].rr) /
                                                      (float) (*ctx_avrQRSRR) >
                                                      0.95) && ((float) (MultiQrstFeatures[peakNumberPre].rr +
                                                                         MultiQrstFeatures[peakLoc].rr) /
                                                                (float) (*ctx_avrQRSRR) < 1.07));
                            if ((qrsIntervalStrong >= 0.9) && (((qrsWidthNarrow >= 0.93) && m) ||
                                                               (MultiQrstFeatures[peakNumberPre].Psum /
                                                                (*ctx_currentPsum) >
                                                                1.8) ||
                                                               ((MultiQrstFeatures[peakNumberPre].prInterval) >= 110) ||
                                                               (MultiQrstFeatures[peakNumberPre].pType == 3) ||
                                                               (MultiQrstFeatures[peakNumberPre].pType == 5)))//%(
                                //%   if((qrsFeatures(peakNumberPre).pType~=1)&&(qrsFeatures(peakNumberPre).pType~=6)||(morphologyIndex>0.98)||(ectopicIndex==1)||(m&&(qrsFeatures(peakNumberPre).pType~=6)))
                            {
                                if (((MultiQrstFeatures[peakNumberPre].pType != 1) &&
                                     (MultiQrstFeatures[peakNumberPre].pType != 6)) || (ectopicIndex == 1) ||
                                    (m && (MultiQrstFeatures[peakNumberPre].pType != 6))) {
                                    if ((MultiQrstFeatures[peakNumberPre].pType == 2) &&
                                        (MultiQrstFeatures[peakNumberPre].prInterval < BEAT_MS120))
                                        (*ctx_japc) = (*ctx_japc) + 1;
                                    else {
                                        (*ctx_apc) = (*ctx_apc) + 1;
                                    }//end
                                    apcFlag = 1;
                                    (*ctx_complementFlag) = 1;
                                }//end
                            } else {
                                if ((MultiQrstFeatures[peakNumberPre].pType != 1) &&
                                    (MultiQrstFeatures[peakNumberPre].qrsWidth < BEAT_MS120)
                                    && (MultiQrstFeatures[peakNumberPre].pType != 6)) {
                                    if ((~interpolatedIndex) && ((MultiQrstFeatures[peakNumberPre].pType == 4) ||
                                                                 (MultiQrstFeatures[peakNumberPre].pType == 3) ||
                                                                 (MultiQrstFeatures[peakNumberPre].pType == 0) ||
                                                                 ((*ctx_currentPsum) /
                                                                  (MultiQrstFeatures[peakNumberPre].Psum) > 1.35))) {
                                        (*ctx_apc) = (*ctx_apc) + 1;
                                        apcFlag = 1;
                                        (*ctx_complementFlag) = 1;
                                    }//end
                                }//end
                            }//end

                            //% sinusParameter=((qrsFeatures(peakNumberPre).Pwidth)<BEAT_MS120)&&(qrsFeatures(peakNumberPre).prInterval<2*BEAT_MS100)&&(pArea>0.85)&&(pArea<1.15);
                            // %  if(((qrsAreaRation<1.16))&&(qrsIntervalStrong==1)&&(qrsFeatures(peakNumberPre).prInterval>2*BEAT_MS110)&&((qrsFeatures(peakNumberPre).qrsWidth)<BEAT_MS110)&&(apcFlag==0))  % BEAT_MS110   0.75 qrsIntervalStrong>=0.75
                            int rBBBLogic = ((morphologyIndex > 0.98) ||
                                             ((MultiQrstFeatures[peakNumberPre].prInterval) >= 2 * BEAT_MS110)) &&
                                            ((MultiQrstFeatures[peakNumberPre].qrsWidth) >= BEAT_MS140) &&
                                            ((qrsAreaRation < 1.26) ||
                                             ((MultiQrstFeatures[peakNumberPre].prInterval) >= 2 * BEAT_MS110)) &&
                                            (apcFlag == 0);
                            if (((qrsAreaRation < 1.16)) && (qrsIntervalStrong == 1) &&
                                ((MultiQrstFeatures[peakNumberPre].qrsWidth) < BEAT_MS110) &&
                                (apcFlag == 0))//%||rBBBLogic==1)
                            {
                                (*ctx_apc) = (*ctx_apc) + 1;
                                apcFlag = 1;
                                (*ctx_complementFlag) = 1;
                            }//end
                        }// end      %%   if((pvcFlag==0)&&(morphologyIndex>=0.92)&&(qrsIntervalStrong==1)) %%%%%%%%%%%%  APC %%%%%%%%%%%%


                        //  %%% 11/12 2022
                    }//   end%% 2512  if(((qrsIntervalStrong>0.749)||(ventricularIndex==1)&&(qrsParameter==1))&&(japcFlag==0))


                    // %%%%%%%%%%%%%%%%%%%% Search back %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ///%  if((peakNumberPre>100)&&( morphologyIndex<0.90)&&(qrsParameter==1)&&(ll<(*ctx_currentWidthMax))&&tWaveDirection)
                    //% if(( morphologyIndex<0.90)&&(pvcAfFlag==1)&&(qrsParameter==1)&&(pvcFlag==0)&&(apcFlag==0))
                    if ((morphologyIndex < 0.90) && (qrsParameter == 1) && (pvcFlag == 0) && (apcFlag == 0) &&
                        (MultiQrstFeatures[peakNumberPre].qrsWidth) > BEAT_MS140) {
                        int pvcExceptEx = ((((MultiQrstFeatures[peakNumberPre].tAera) / (*ctx_currentTArea)) > 2.5) &&
                                           (((MultiQrstFeatures[peakNumberPre].speak) / (*ctx_currentAverSpeak)) > 8));
                        int pvcExtraLogic = ((((MultiQrstFeatures[peakNumberPre].pPeak) < 3.0) &&
                                              (((MultiQrstFeatures[peakNumberPre].tAera) / (*ctx_currentTArea)) > 2.5 ||
                                               (MultiQrstFeatures[peakNumberPre].qrsWidth) > BEAT_MS140 ||
                                               tWaveDirection ==
                                               1)));//%%||(qrsFeatures(peakNumberPre).qrsWidth)>BEAT_MS160);
                        if (((fuzzyNormalResult <= nomalRation)) &&
                            (((qrsWidthWide == 1) && (qrsIntervalStrong == 1) && (qrsComplementStrong < pvcThrehold)) ||
                             ((qrsWidthWide == 1) && (qrsIntervalStrong < 0.5) && (qrsComplementStrong == 1))\
 || pvcExceptEx || pvcExtraLogic || ((qrsWidthWide < 0.5) && (qrsIntervalStrong == 1) && (qrsComplementStrong == 1)) ||
                             (MultiQrstFeatures[peakNumberPre].area > 2.5 * (*ctx_currentArea)) ||
                             (((fabsf(MultiQrstFeatures[peakNumberPre].area)) > 1.30 * (*ctx_currentArea)) &&
                              ((qrsIntervalStrong > 0.93)))))//%%(((abs(qrsFeatures((*ctx_peakNumber)-1).speak
                        {
                            //%    ||((qrsWidthWide<0.5)&&(qrsIntervalStrong==1)&&(qrsComplementStrong==1)) ||(qrsFeatures((*ctx_peakNumber)-1).area>2.5*ctx_currentArea)||(((abs(qrsFeatures((*ctx_peakNumber)-1).area))>1.30*ctx_currentArea)&&((qrsIntervalStrong>0.93)))))%%(((abs(qrsFeatures((*ctx_peakNumber)-1).speak
//%   ||((qrsWidthWide<0.5)&&(qrsIntervalStrong==1)&&(qrsComplementStrong==1))||(qrsFeatures((*ctx_peakNumber)-1).qrsWidth>=75) ||(qrsFeatures((*ctx_peakNumber)-1).area>2.5*ctx_currentArea)||(((abs(qrsFeatures((*ctx_peakNumber)-1).area))>1.30*ctx_currentArea)&&((qrsIntervalStrong>0.93)))))%%(((abs(qrsFeatures((*ctx_peakNumber)-1).speak
                            //%  if(qrsIntervalStrong<0.2)&&(ventricularIndex==1)&&((abs(qrsFeatures(peakNumberPre).pPeak))>(abs(qrsFeatures(peakNumberPre).qrsMainWave))/10)%% for p wave exception
                            if ((qrsIntervalStrong < 0.2) &&
                                ((fabsf(MultiQrstFeatures[peakNumberPre].pPeak)) > (*ctx_currentPeak) * 0.8)) {
                                fuzzyPVCResult = 0;
                            }//end
                            int m;
                            if ((fuzzyPVCResult > 0) && (apcFlag == 0)) {
                                m = 0;
                                while (m < 100) {
                                    newArray[m] = lastBeat[m + 170];
                                    m = m + 1;
                                }//end
                                float maxCoef = templateSeek(&(*ctx_pvcLibInfo), newArray);
                                if (maxCoef >
                                    0.92) {//%%0.92))%%&&(qrsFeatures(peakNumberPre).area)<12000)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test15:09 06/12 20181206 ypz   %0.92)
                                    (*ctx_pvcNumber) = (*ctx_pvcNumber) + 1;
                                    (*ctx_complementFlag) = 1;
                                    pvcFlag = 1;
                                }//end

                            }//end
                        }//end
                    }//end
//%%%%%%%   20221203  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//% end %% end if >=8    Line 1074


// %%%%%%%%%%%% 09:38 13/12 2022  %%%%%%%%%%%%%%%%%%%%%%%%%%

                    if (((MultiQrstFeatures[peakNumberPre].qrsWidth) > 65) && (apcFlag == 0) && (pvcFlag == 0)) {

                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%% LBBB RBBB  %%%%%%%%%%%%%%%%%%%%%%%%%%BEAT_MS140

                        if (MultiQrstFeatures[peakNumberPre].morphologyType == 2) {
                            (*ctx_rbbbCount) = (*ctx_rbbbCount) + 1;
                        } else {
                            if (MultiQrstFeatures[peakNumberPre].morphologyType == 1) {
                                (*ctx_lbbbCount) = (*ctx_lbbbCount) + 1;
                            }//end
                        }//end
                        //%%%%%%%%%%%%%%%  WPW  %%%%%%%%%%%%%%%%%%%%%%;
                        if ((MultiQrstFeatures[peakNumberPre].qPeak == 0) &&
                            (MultiQrstFeatures[peakNumberPre].notchFlag == 0)) {
                            float temp = MultiQrstFeatures[peakNumberPre].qrsHeight;
                            int pwpResult = wpwCalculate( lastBeat, MultiQrstFeatures[peakNumberPre].onset, temp);
                            if (pwpResult == 1) {
                                (*ctx_wpwCount) = (*ctx_wpwCount) + 1;

                            }//end
                        }//end
                    }//end

                    //%%%%%%%%%%%%%%%%%%%%%%%%%%
                    //%apcCount( apcFlag ,0);
                    //  %%%%%%   ventricular bigeminy,trigeminy,tachycardia    %%%%%%%%%%% apcFlag
                    ventricularPrematureCount( ventricularPrematureCountCtx, pvcFlag,
                                              MultiQrstFeatures[peakNumberPre].rr, &(*ctx_doublePvc),
                                              &(*ctx_vpcBigeminyCount), &(*ctx_vpcTrigeminyCount),
                                              &(*ctx_ventricularTachycardiaCount),
                                              &(*ctx_longestPvcTachycardia), &(*ctx_fastestPvcTachy), 0);
                    // %%%%%%%%%%%%%%%%%%%%%%%%%%
                    atrialPrematureCount( atrialPrematureCountCtx, apcFlag, MultiQrstFeatures[peakNumberPre].rr,
                                         &(*ctx_doubleAvc), &(*ctx_apcBigeminyCount),
                                         &(*ctx_apcTrigeminyCount), &(*ctx_atrialTachycardiaConut),
                                         &(*ctx_longestApcTachycardia),
                                         &(*ctx_fastestApcTachy), 0);

                    //%%%%%%%%%%%%%%%%%%%%%% escape  %%%%%%%%%%%%%%%%%%%%%%%%%%
                    if ((MultiQrstFeatures[peakNumberPre].pType == 2) && (MultiQrstFeatures[peakNumberPre].rr >= 630) &&
                        (MultiQrstFeatures[peakNumberPre].prInterval < BEAT_MS150) &&
                        (morphologyIndex >= 0.91))//%&&((*ctx_currentPeak)/abs(qrsFeatures(peakNumberPre).pPeak))>1.5)
                        // junctionalEscapeCount=junctionalEscapeCount+1;
                        //end
                        //%%%%%%%%%%%%%%%%%%  ST Segment algorthm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if ((((*ctx_complementFlag) == 0) || ((*ctx_complementFlag) == 3)) &&
                            ((*ctx_learningNumber) >= templateLength) &&
                            (((*ctx_learningFlag) == 1) && (japcFlag == 0))) {
                            //[ ST_MorphologyType, ST_Change,ST_ElavatorAmph,type] = ST_SegmentAnalysis( lastBeat, MultiQrstFeatures(peakNumberPre).rr, floor(MultiQrstFeatures(peakNumberPre).Poffset) ,floor(MultiQrstFeatures(peakNumberPre).Ponset) ,floor(MultiQrstFeatures(peakNumberPre).onset),floor(MultiQrstFeatures(peakNumberPre).offset),amplifierGainConst,MultiQrstFeatures(peakNumberPre).tFirstPeak,Crit,leadStandard);
                            //[ ST_MorphologyType, ST_Change,ST_ElavatorAmph,type] = ST_SegmentAnalysis( lastBeat, qrsFeatures(peakNumberPre).rr, qrsFeatures(peakNumberPre).Ponset ,qrsFeatures(peakNumberPre).Poffset ,qrsFeatures(peakNumberPre).onset,qrsFeatures(peakNumberPre).offset,amplifierGainConst,qrsFeatures(peakNumberPre).tFirstPeak,Crit,leadStandard);
                            ST_SegmentAnalysis( lastBeat, MultiQrstFeatures[peakNumberPre].rr,
                                               MultiQrstFeatures[peakNumberPre].Ponset,
                                               MultiQrstFeatures[peakNumberPre].Poffset,
                                               MultiQrstFeatures[peakNumberPre].onset,
                                               MultiQrstFeatures[peakNumberPre].offset, amplifierGainConst,
                                               MultiQrstFeatures[peakNumberPre].tFirstPeak, Crit, leadStandard,
                                               &ST_MorphologyType, &ST_Change, &ST_ElavatorAmph, &type);

                            MultiQrstFeatures[peakNumberPre].stChange = (int) ST_ElavatorAmph;
                            MultiQrstFeatures[peakNumberPre].stType = type;
                            int cycle = 1;
                            int ST_ChangeCount = 0;
                            int ST_ChangeFlag = 0;
                            if (fabsf(ST_ElavatorAmph) > 0.12) {
                                (*ctx_stChangeCount) = (*ctx_stChangeCount) + 1;
                            }//end
                            while (cycle <= 20) {
                                ST_IschaemicEpisode[cycle] = ST_IschaemicEpisode[cycle + 1];
                                if (ST_IschaemicEpisode[cycle] == 1) {
                                    ST_ChangeCount = ST_ChangeCount + 1;
                                }//end
                                cycle = cycle + 1;
                            }//end

                            if (ST_Change != 0) {
                                ST_IschaemicEpisode[19] = 1;
                                (*ctx_stChangeCount) = (*ctx_stChangeCount) + 1;
                            }//end

                            if (((float) ST_ChangeCount / 20) > 0.75) {
                                ST_ChangeFlag = 1;
                            }//end

                        }//end
                    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normal QRS parameters update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    int peakTempBufferFront = ((*ctx_MultiPeakNumber) -
                                               2);                           //%%%%%%%%%%%%%%%%%%%
                    if (peakTempBufferFront < 0)
                        peakTempBuffer = peakTempBuffer + maxQueueLength;
                    //end
                    int peakTempBuffersecond = ((*ctx_MultiPeakNumber) -
                                                1);                           //%%%%%%%%%%%%%%%%%%%
                    if (peakTempBuffersecond < 1)
                        peakTempBuffer = peakTempBuffer + maxQueueLength;
                    //end

                    int lbbbHRVLogic, lbbbWidthUpdateLogic, lbbbUpdateLogic;
                    if ((*ctx_MultiPeakNumber) > 5) {
                        lbbbHRVLogic = ((abs(MultiQrstFeatures[peakTempBufferFront].rr) -
                                         MultiQrstFeatures[peakTempBuffersecond].rr) < BEAT_MS140 / 2) &&
                                       (((abs(MultiQrstFeatures[peakTempBuffersecond].rr) -
                                          MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr)) < BEAT_MS140 / 2);
                        lbbbWidthUpdateLogic = ((MultiQrstFeatures[peakTempBufferFront].rr > BEAT_MS500) &&
                                                (MultiQrstFeatures[peakTempBuffersecond].rr > BEAT_MS500) &&
                                                (MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr > BEAT_MS500));
                        lbbbUpdateLogic = ((MultiQrstFeatures[peakTempBufferFront].qrsWidth > BEAT_MS140) &&
                                           (MultiQrstFeatures[peakTempBuffersecond].qrsWidth > BEAT_MS140) &&
                                           (MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsWidth > BEAT_MS140));
                    } else {
                        lbbbHRVLogic = 0;
                        lbbbWidthUpdateLogic = 0;
                        lbbbUpdateLogic = 0;
                    }//end
                    int LBBBUpdateLogic = lbbbUpdateLogic && lbbbWidthUpdateLogic && lbbbHRVLogic;

                    //%%%%%%%%%%%% sinus or lbbb or rbbb update logic %%%%%%%%%%
                    int kk = 0;
                    if (((fuzzyNormalResult > 0.8) && (areaIndex < 1.35) && (japcFlag == 0) &&
                         ((*ctx_learningNumber) >= templateLength)) ||
                        ((fuzzyNormalResult == 0) && (LBBBUpdateLogic == 1))) {
                        kk = 0;
                        while (kk < 3) {
                            qrsSpeak[kk] = qrsSpeak[kk + 1];
                            qrsArea[kk] = qrsArea[kk + 1];
                            tWaveArea[kk] = tWaveArea[kk + 1];

                            peakTempBuffer = ((*ctx_MultiPeakNumber) - 1);

                            if (peakTempBuffer < 0) peakTempBuffer = peakTempBuffer + maxQueueLength;

                            if ((MultiQrstFeatures[peakTempBuffer].pPeak) != 0) {
                                Pheight[kk] = Pheight[kk + 1];
                                Parea[kk] = Parea[kk + 1];
                                tNormalAverage[kk] = tNormalAverage[kk + 1];
                            }//end
                            kk = kk + 1;
                        }//end
                        peakTempBuffer = ((*ctx_MultiPeakNumber) - 1);
                        if (peakTempBuffer < 0) peakTempBuffer = peakTempBuffer + maxQueueLength;
                        //end

                        qrsSpeak[kk] = fabsf(MultiQrstFeatures[peakTempBuffer].speak);
                        qrsSpeak[kk] = FeaturesOptimize( (*ctx_currentAverSpeak),
                                                        fabsf(MultiQrstFeatures[peakNumberPre].speak));
                        qrsArea[kk] = FeaturesOptimize( (*ctx_currentArea),
                                                       fabsf(MultiQrstFeatures[peakNumberPre].area));
                        tWaveArea[kk] = FeaturesOptimize((*ctx_currentTArea),
                                                         fabsf((MultiQrstFeatures[peakNumberPre].tAera)));
                        tNormalAverage[kk] = FeaturesOptimize((*ctx_tHeightAverage),
                                                              fabsf((MultiQrstFeatures[peakNumberPre].tHeight)));

                        if ((MultiQrstFeatures[peakNumberPre].pPeak) != 0) {
                            Pheight[kk] = fabsf(MultiQrstFeatures[peakNumberPre].pPeak);
                            Parea[kk] = fabsf(MultiQrstFeatures[peakNumberPre].Psum);
                        }//end


                        (*ctx_currentTArea) = templateCompute(tWaveArea, templateLength);
                        (*ctx_currentAverSpeak) = templateCompute(qrsSpeak, templateLength);
                        (*ctx_currentArea) = templateCompute(qrsArea, templateLength);
                        (*ctx_currentPeak) = templateCompute(Pheight, templateLength);
                        (*ctx_currentPsum) = templateCompute(Parea, templateLength);
                        (*ctx_tHeightAverage) = templateCompute(tNormalAverage, templateLength);
                        //normalTemplate=[];
                        if ((MultiQrstFeatures[peakNumberPre].onset != 0) &&
                            (MultiQrstFeatures[peakNumberPre].offset != 0)) {
                            //normalTemplate= templateCopy(lastBeat ,MultiQrstFeatures[peakNumberPre].onset,MultiQrstFeatures[peakNumberPre].offset);
                            //(*ctx_templateFlag)=1;
                            templateCopy(lastBeat, normalTemplate);
                            (*ctx_templateFlag) = 1;
                        }//end
                        // latestPwave=normalP;
                        // normalP=[];
                        //% normalP=arrayCopy(lastBeat ,qrsFeatures((*ctx_peakNumber)-1).Ponset,qrsFeatures((*ctx_peakNumber)-1).Poffset);

                    }//end %  (*ctx_currentPsum)

                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    // rrIntervalUpdateLogic=0;
                    int rrIntervalUpdateLogic = 0;
                    int rrIntervalTemp;
                    float coef = 0;
                    int weight = 0;
                    float averRR = 0;
                    if (((*ctx_avrQRSRR) == 0) && ((*ctx_learningFlag) == 0)) {
                        (*ctx_avrQRSRR) = (MultiQrstFeatures[peakNumberPre].rr);
                    }//end

                    if (apcFlag == 1 || pvcFlag == 1)  //%% complementary beats
                    {
                        rrIntervalTemp =
                                ((MultiQrstFeatures[peakNumberPre].rr) +
                                 (MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr)) / 2;
                        coef = (float) rrIntervalTemp / (float) (*ctx_avrQRSRR);

                        if ((coef > 0.91) && (coef < 1.09)) {
                            rrIntervalUpdateLogic = 1;
                        }//end
                    }// end

                    if ((((*ctx_complementFlag) == 0) || ((*ctx_complementFlag) == 3) || rrIntervalUpdateLogic == 1) &&
                        (((*ctx_learningFlag) == 1) && (japcFlag == 0))) {

                        if (rrIntervalUpdateLogic == 0) {
                            coef = (float) MultiQrstFeatures[peakNumberPre].rr / (float) (*ctx_avrQRSRR);
                        }// end
                        //%    coef=(qrsFeatures(peakNumberPre).rr)/(*ctx_avrQRSRR);
                        //int k;
                        if ((coef >= 0.5) && (coef <= 1.5)) {  //%% only within [50% 180%] can be updated
                            kk = 0;
                            while (kk < (templateLength - 1)) {
                                qrsRR[kk] = qrsRR[kk + 1];
                                qrsWidth[kk] = qrsWidth[kk + 1];
                                kk = kk + 1;
                            }

                            if ((coef >= 0.85) && (coef <= 1.15)) {
                                weight = 0;
                            } else {
                                if (((coef >= 0.7) && (coef <= 0.9)) || ((coef >= 1.1) && (coef <= 1.5))) {
                                    weight = 3;//%3;
                                } else {
                                    if (((coef >= 0.4) && (coef <= 0.7)) || ((coef >= 1.5) && (coef <= 1.8))) {
                                        weight = 7;
                                    } else {
                                        weight = 8;
                                    }//end
                                }//end
                            }// end


                            if (rrIntervalUpdateLogic == 1)
                                averRR = (float) rrIntervalTemp;
                            else {
                                averRR = (float) MultiQrstFeatures[peakNumberPre].rr;
                            }// end

                            if ((coef >= 0.85) && coef <= 1.15) {
                                qrsRR[kk] = averRR;
                            } else {
                                // %     qrsRR(k)=((*ctx_avrQRSRR)+averRR*weight)/(weight+1);
                                //%        qrsRR(k)=((*ctx_avrQRSRR)*weight+averRR)/(weight+1);
                                if (weight <= 3) {
                                    qrsRR[kk] =
                                            ((float) (*ctx_avrQRSRR) + averRR * (float) weight) / (float) (weight + 1);
                                } else {
                                    qrsRR[kk] = ((float) ((*ctx_avrQRSRR) * weight) + averRR) / (float) (weight + 1);
                                }//end
                            }//end
                            (*ctx_avrQRSRR) = (int) templateCompute(qrsRR, templateLength);
                        }
                    }//end if((((*ctx_complementFlag)==0)||((*ctx_complementFlag)==3)|| rrIntervalUpdateLogic==1)&&(((*ctx_learningFlag)==1)&&(japcFlag==0)))


                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   14:29 11/05 2019  ypz

                    result->sinus_arrhythmia.sinusArrestType = (*ctx_sinusArrestType);
                    result->sinus_arrhythmia.sinusArrhythmiaType = (*ctx_sinusArrhythmiaType);
                    result->sinus_arrhythmia.sinusBradycardiaType = (*ctx_sinusBradycardiaType);
                    result->sinus_arrhythmia.sinusTachycardiaType = (*ctx_sinusTachycardiaType);
                    result->sinus_arrhythmia.sickSinusSyndromeType = (*ctx_sickSinusSyndromeType);
                    result->sinus_arrhythmia.sinusIrregularityType = (*ctx_sinusIrregularityType);

                    if (apcFlag == 1) {
                        MultiQrstFeatures[peakNumberPre].type = 1;
                        (*ctx_apcAFlag) = 1;
                        result->supra_arrhythmia.apcTotal = (*ctx_apc);
                        result->supra_arrhythmia.doubleApc = (*ctx_doubleAvc);
                        result->supra_arrhythmia.fastestApcTachy = (*ctx_fastestApcTachy);
                        result->supra_arrhythmia.apcBigeminyCount = (*ctx_apcBigeminyCount);
                        result->supra_arrhythmia.longestApcTachycardia = (*ctx_longestApcTachycardia);
                        result->supra_arrhythmia.atrialTachycardiaConut = (*ctx_atrialTachycardiaConut);
                        result->supra_arrhythmia.apcTrigeminyCount = (*ctx_apcTrigeminyCount);
                    } else {
                        if (pvcFlag == 1) {
                            MultiQrstFeatures[peakNumberPre].type = 2;
                            //(*ctx_pvcAFlag)=1;
                            //(*ctx_pvcLength)= (*ctx_pvcLength)+1;
                            //pvcIndex=1;
                            //qrsAverageWidth =0;
                            //pvcAverage=0;
                            (*ctx_pvcAFlag) = 1;
                            (*ctx_pvcLength) = (*ctx_pvcLength) + 1;
                            int pvcIndex = 0;
                            int qrsAverageWidth = 0;
                            float pvcAverage = 0;
                            float ll;
                            if (((*ctx_pvcLength) <= pvcQueLength) && (MultiQrstFeatures[peakNumberPre].area != 0)) {
                                // % qrsFeatures(peakNumberPre).qrsWidth;
                                while (pvcIndex < pvcQueLength) {
                                    pvcWidthArray[pvcIndex] = pvcWidthArray[pvcIndex + 1]; //% pvcQrsWidthArray
                                    pvcQrsWidthArray[pvcIndex] = pvcQrsWidthArray[pvcIndex + 1];
                                    pvcIndex = pvcIndex + 1;
                                }//end
                                pvcWidthArray[pvcQueLength - 1] = (int) MultiQrstFeatures[peakNumberPre].area; //%
                                pvcQrsWidthArray[pvcQueLength - 1] = MultiQrstFeatures[peakNumberPre].qrsWidth; //%
                                ll = MultiQrstFeatures[peakNumberPre].area;

                            }//end
                            result->ventricular_arrhythmia.pvcCount = (*ctx_pvc) + (*ctx_pvcNumber);
                            result->ventricular_arrhythmia.doublePvc = (*ctx_doublePvc);
                            result->ventricular_arrhythmia.fastestPvcTachy = (*ctx_fastestPvcTachy);
                            result->ventricular_arrhythmia.vpcBigeminyCount = (*ctx_vpcBigeminyCount);
                            result->ventricular_arrhythmia.vpcTrigeminyCount = (*ctx_vpcTrigeminyCount);
                            result->ventricular_arrhythmia.longestPvcTachycardia = (*ctx_longestPvcTachycardia);
                            result->ventricular_arrhythmia.ventricularTachycardiaCount = (*ctx_ventricularTachycardiaCount);

                            if (((*ctx_pvcLength) > pvcQueLength) &&
                                ((MultiQrstFeatures[peakNumberPre].area) < (float) (*ctx_currentWidthMax)) &&
                                (MultiQrstFeatures[peakNumberPre].qrsWidth < (*ctx_currentQRS)) &&
                                (MultiQrstFeatures[peakNumberPre].area != 0)) {
                                while (pvcIndex < pvcQueLength) {
                                    pvcWidthArray[pvcIndex] = pvcWidthArray[pvcIndex + 1];
                                    pvcQrsWidthArray[pvcIndex] = pvcQrsWidthArray[pvcIndex + 1];
                                    pvcIndex = pvcIndex + 1;
                                }//end
                                pvcWidthArray[pvcQueLength - 1] = (int) MultiQrstFeatures[peakNumberPre].area;
                                pvcQrsWidthArray[pvcQueLength - 1] = MultiQrstFeatures[peakNumberPre].qrsWidth;
                                ll = MultiQrstFeatures[peakNumberPre].area;
                            }//end

                            if ((*ctx_pvcLength) >= pvcQueLength) {
                                pvcIndex = 1;
                                while (pvcIndex <= pvcQueLength) {
                                    qrsAverageWidth = qrsAverageWidth + pvcWidthArray[pvcIndex];
                                    pvcAverage = pvcAverage + (float) pvcQrsWidthArray[pvcIndex];
                                    pvcIndex = pvcIndex + 1;
                                }//end
                                float averageArea = (float) qrsAverageWidth / pvcQueLength;
                                float diffArea = 0;
                                pvcIndex = 1;
                                while (pvcIndex <= pvcQueLength) {
                                    diffArea = diffArea + ((float) pvcWidthArray[pvcIndex] - averageArea) *
                                                          ((float) pvcWidthArray[pvcIndex] - averageArea);
                                    pvcIndex = pvcIndex + 1;
                                }//end
                                float llll = sqrtf(diffArea / pvcQueLength) / averageArea;
                                (*ctx_currentQRS) = (int) (pvcAverage * 1.3 / pvcQueLength);
                                if (llll < 0.15) {
                                    (*ctx_currentWidthMax) = (int) ((double) qrsAverageWidth * 1.55 /
                                                                    pvcQueLength); //% 1.5/6;
                                } else {
                                    (*ctx_currentWidthMax) = (int) ((double) qrsAverageWidth * 3.0 / pvcQueLength);
                                }//end
                            }// end if((*ctx_pvcLength)>=pvcQueLength)

                        }//end    if(pvcFlag==1)
                    }//end  else  if



                    //%%%%%%%%%%%%%%%%%% atria and ventricular   escape  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    int junctionalLogics;
                    int tPvcDirectory = ((((MultiQrstFeatures[peakNumberPre].qrsMainWave) *
                                           (MultiQrstFeatures[peakNumberPre].tHeight)) <
                                          0));//%&&(morphologyIndex<0.91));
                    int ventricularEscape = ((fabsf(MultiQrstFeatures[peakNumberPre].pPeak) < (*ctx_currentPeak) / 3) &&
                                             ((MultiQrstFeatures[peakNumberPre].qrsWidth) >= BEAT_MS140) &&
                                             (tPvcDirectory) &&
                                             (MultiQrstFeatures[peakNumberPre].area > 2.0 * (*ctx_currentArea)));
                    if (((*ctx_complementFlag) != 2) && (MultiQrstFeatures[peakNumberPre].rr >= BEAT_MS1100) &&
                        (MultiQrstFeatures[peakNumberPre].rr <= BEAT_MS3000) && (ventricularEscape == 1)) {
                        (*ctx_ventrEscapeCount) = (*ctx_ventrEscapeCount) + 1;
                    } else
                        junctionalLogics = (((fabsf(MultiQrstFeatures[peakNumberPre].pPeak) < (*ctx_currentPeak) / 3) &&
                                             ((MultiQrstFeatures[peakNumberPre].prInterval) <= BEAT_MS140)) ||
                                            (((MultiQrstFeatures[peakNumberPre].prInterval) <= BEAT_MS120) &&
                                             (MultiQrstFeatures[peakNumberPre].pType) == 2));
                    int junctionalEscape = ((morphologyIndex > 0.93) &&
                                            ((MultiQrstFeatures[peakNumberPre].qrsWidth) <= BEAT_MS120) &&
                                            ((MultiQrstFeatures[peakNumberPre].prInterval) <= BEAT_MS120 ||
                                             junctionalLogics));
                    if (((*ctx_complementFlag) == 0) && (MultiQrstFeatures[peakNumberPre].rr >= BEAT_MS1000) &&
                        (MultiQrstFeatures[peakNumberPre].rr <= BEAT_MS4000) && (junctionalEscape == 1))
                        (*ctx_junctEscapeCount) = (*ctx_junctEscapeCount) + 1;
                    else {
                        if ((MultiQrstFeatures[peakNumberPre].prInterval) > BEAT_MS120) {
                            //%  ectopicEscapeCount=ectopicEscapeCount+1;
                        }// end
                    }//end
                    //end




//%%%%%%%%%%%%%%%%%%%%%%%% Sinus AV and intraventricular conduction block   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    // [ AVBlockFlag] = firstDegreeAVBlock((MultiQrstFeatures(peakNumberPre).prInterval),0);
                    int AVBlockFlag = firstDegreeAVBlock( firstDegreeAvBlockCtx,
                                                         MultiQrstFeatures[peakNumberPre].prInterval, 0);
//%   [sinusPwaveCount] = pWaveFind(array,rr,(*ctx_currentPeak));
                    if ((MultiQrstFeatures[peakNumberPre].ppInterval != 0)) {
                        if ((*ctx_complementFlag) != 0) {
                            // sinoAtrialBlock(0,0);
                        } else {
                            sinusArrhythmia( sinusArrhythmiaCtx, (MultiQrstFeatures[peakNumberPre].ppInterval),
                                            &(*ctx_sinusArrhythmiaType),
                                            &(*ctx_sinusIrregularityType), &(*ctx_sinusArrestType),
                                            &(*ctx_sinusBradycardiaType),
                                            &(*ctx_sinusTachycardiaType), &(*ctx_sickSinusSyndromeType), 0);
                        }
                        //%  [sinoAtrialBlockType] = sinusArrhythmia((qrsFeatures(peakNumberPre).ppInterval),0);
                        //[(*ctx_sinusArrhythmiaType),(*ctx_sinusIrregularityType),(*ctx_sinusArrestType),~,(*ctx_sinusTachycardiaType),sickSinusSyndromeType] = sinusArrhythmia((floor(MultiQrstFeatures(peakNumberPre).ppInterval)),0);
//%         sinoAtrialBlock((qrsFeatures(peakNumberPre).ppInterval),0);
                        //end
                    } //  end  if((MultiQrstFeatures[peakNumberPre].ppInterval!=0))
//end

//%%%%%%%%%%%%%%%%%%%%%%%%%%% long rr interval    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    int jjEnd, jjStart;
                    int rrPeriod = MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr;

                    if ((rrPeriod > BEAT_MS1500) && (rrPeriod < 2 * BEAT_MS1000)) {
                        jjEnd = Interval;
                        jjStart = Interval - MultiQrstFeatures[(*ctx_MultiPeakNumber)].rr;
                        if (jjStart < 1) {
                            jjStart = jjStart + peakLocationMod;
                        }//end

                        int arrayStart = jjStart + (MultiQrstFeatures[peakNumberPre].tOffset - FIDMARK);
                        if (arrayStart > peakLocationMod) {
                            arrayStart = arrayStart - peakLocationMod;
                        }//end

                        int arrayEnd = jjEnd - (FIDMARK - MultiQrstFeatures[(*ctx_MultiPeakNumber)].onset);
                        if (arrayEnd < 1) {
                            arrayEnd = arrayEnd + peakLocationMod;
                        }//end

                        len = arrayEnd - arrayStart;
                        if (len < 1) {
                            len = len + peakLocationMod;
                        }//end

                        int jj = arrayStart;
                        int kk = 0;
                        float tempArray[len];//=zeros(1,len);

                        while (kk < len) {
                            tempArray[kk] = ECGBuffer[2][jj];
                            kk = kk + 1;
                            jj = jj + 1;
                            if (jj > peakLocationMod) {
                                jj = jj - peakLocationMod;
                            }//end
                        }//end


                        start = Interval - (FIDMARK + (int) ceil(0.06 * sampleRate));
                        if (start < 1) {
                            start = start + peakLocationMod;
                        }//end
                        finish = Interval + (int) ceil(0.06 * sampleRate);
                        if (finish > peakLocationMod) {
                            finish = finish - peakLocationMod;
                        }//end
                        //%      [ AAA ] = bufferCopy( ECGBuffer ,start,finish);
                        //%
                        //%       [index] = checkQuality(AAA)
                        //%       if(index<0)

                        // int peakNumber1= longRRBlock(tempArray,len,MultiQrstFeatures[(*ctx_MultiPeakNumber)].qrsHeight,(*ctx_currentPeak),(*ctx_currentPsum));
                        // MultiQrstFeatures[(*ctx_MultiPeakNumber)].blockedNumber=peakNumber1;
                        // if( peakNumber1>1){

                        // (*ctx_twoDgreeBlockNumber)=(*ctx_twoDgreeBlockNumber)+1;
                        // };//end
                    }// end if((rrPeriod>BEAT_MS1500)&&(rrPeriod<2*BEAT_MS1000))


                    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                }//end %% end if >=8 	        Line 4110


            }//end %% if(rrCount>=3)         Line 3692


        }// if(rrValue!=0)            Line 3678
//        else {
//            jj++;
//        }

    }//	 while(jj<ChannelNo)           Line 3668





    return 0;

}//  int ArrhythmiaAnalysis20230311  line 3197




void init_arrhythmia_type(arrhythmia_type *at) {
    at->supra_arrhythmia.apcTotal = 0;
    at->supra_arrhythmia.japcTotal = 0;
    at->supra_arrhythmia.doubleApc = 0;
    at->supra_arrhythmia.apcBigeminyCount = 0;
    at->supra_arrhythmia.apcTrigeminyCount = 0;
    at->supra_arrhythmia.atrialTachycardiaConut = 0;
    at->supra_arrhythmia.longestApcTachycardia = 0;
    at->supra_arrhythmia.fastestApcTachy = 0;
    at->supra_arrhythmia.atrial_fibrillation = 0;
    at->supra_arrhythmia.atrial_flutter = 0;
    at->ventricular_arrhythmia.doublePvc = 0;
    at->ventricular_arrhythmia.fastestPvcTachy = 0;
    at->ventricular_arrhythmia.lbbbCount = 0;
    at->ventricular_arrhythmia.longestPvcTachycardia = 0;
    at->ventricular_arrhythmia.pvcCount = 0;
    at->ventricular_arrhythmia.ventricularTachycardiaCount = 0;
    at->ventricular_arrhythmia.rbbbCount = 0;
    at->ventricular_arrhythmia.vpcBigeminyCount = 0;
    at->ventricular_arrhythmia.vpcTrigeminyCount = 0;
    at->block_Count.AsystolyCount = 0;
    at->block_Count.high_degreeCount = 0;
    at->block_Count.one_degreeCount = 0;
    at->block_Count.three_degreeCount = 0;
    at->block_Count.two_degreeCount = 0;
    at->escape_Count.junctional_escape = 0;
    at->escape_Count.ventricular_escape = 0;
    at->sinus_arrhythmia.sickSinusSyndromeType = 0;
    at->sinus_arrhythmia.sinusArrestType = 0;
    at->sinus_arrhythmia.sinusArrhythmiaType = 0;
    at->sinus_arrhythmia.sinusBradycardiaType = 0;
    at->sinus_arrhythmia.sinusIrregularityType = 0;
    at->sinus_arrhythmia.sinusTachycardiaType = 0;
    at->st_value.st_value = 0;
    at->st_value.type = 0;
}

							
			