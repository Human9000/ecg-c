#include "beatFeatureGet.h"
#include "utility.h"
#include <math.h>
#include <string.h>
//#include <QtCore>

#define   FIDMARK        200
#define   BEAT_MS10      5
#define   BEAT_MS8       4
#define   BEAT_MS20      10
#define   BEAT_MS90      45
#define   BEAT_MS250     125
#define   BEAT_MS40      20
#define   arrayLimit     300
#define   BEAT_MS500     250
#define   BEAT_MS6       3
#define   BEAT_MS1000    500
#define   ISO_LIMT       20
#define   EPSILON        0.000001


#define speakThrehold       228
#define Fs                  500
#define dataLength          500
#define BEAT_SAMPLE_RATE    500
#define BEAT_MS_PER_SAMPLE  (1000/BEAT_SAMPLE_RATE)
#define BEAT_MS12           round(12/BEAT_MS_PER_SAMPLE + 0.5)
#define BEAT_MS30           15
#define BEAT_MS50           round (50/BEAT_MS_PER_SAMPLE + 0.5)
#define BEAT_MS80           round (80/BEAT_MS_PER_SAMPLE + 0.5)
#define BEAT_MS400          200
#define BEAT_MS300          round (300/BEAT_MS_PER_SAMPLE + 0.5)
#define BEAT_MS500          250
#define BEAT_MS480          240
#define BEAT_MS1000         500
#define BEATLGTH            BEAT_MS1000
#define BEAT_MS160          80
#define ISO_LENGTH          15
#define BEAT_MS520          260
#define BEAT_MS120          60
#define ISO_LENGTH_DIFF     15
#define BEAT_MS20           10
#define BEAT_MS24           12
#define BEAT_MS10           5
#define BEAT_MS6            3
#define BEAT_MS16           8
#define slopeLength         8
#define TRUE                1
#define FALSE               0
#define slopeLength         8






// qrsFeaturesClass::qrsFeaturesClass()
//  {


//     int kk=0;
//     while(kk<arrayLimit)
//     {
//        firstDiff[kk] =0;
//        secondDiff[kk]=0;
//        fnew[kk]=0;
//        result[kk]=0;
//        data[kk]=0;
//        y[kk]=0;

//        kk++;
//     }

//     kk=0;
//    while(kk<BEAT_MS20)
//    {
//       QRSPeaks[kk]=0;
//       QRSPeaksPositions[kk]=0;
//      kk++;
//    }

//    kk=0;
//   while(kk<sampleRate)
//   {
//      Xpc[kk] =0;
//      kk++;
//   }



//   //////////////////////////////////////////////////////////////////////
// / \brief qrsFeatures::pIsoLevel
// / \param array
// / \param position
// / \param


float pIsoLevelB(const float *array, int position) {
    float averResult = 0;
    float pSum;
    if ((position - 8) > 0) {
        pSum = 0;
        int index = position - 8;
        while (index <= position) {
            pSum = pSum + array[index];
            index = index + 1;
        }//end
        averResult = pSum / 8;
    }//end

    return averResult;

}


//   //////////////////////////////////////////////////////////////////////
// / \brief qrsFeatures::LowPassFilter
// / \param X
// / \param Xpc
// / \param
void LowPassFilterB(const float *X, int fs, int Fpa, float *Xpc) {
    // %  Linear integer lowpass filter

    int mpa = 10;//(int)(fs/Fpa);
    // % c_mpa=20;
    //        sampleRate=500;
    // dataLength=length(X);
//         Xpc=zeros(1,dataLength);
//         data=zeros(1,2*c_mpa);
    float data1[20];
    int kk = 0;
    while (kk < 20)//sampleRate)
    {
        data1[kk] = 0;
        kk++;
    }

    float y0 = 0;
    float y1 = 0;
    float y2 = 0;
    int ptr = 0;
    float output;
    int jj;

    while (ptr < dataLength) {
        y0 = 2 * y1 - y2 + X[ptr] - 2 * data1[mpa - 1] + data1[mpa * 2 - 1];
        y2 = y1;
        y1 = y0;
        output = y0 / (float) (mpa * mpa);
        Xpc[ptr] = output;

        jj = 2 * mpa - 2;
        while (jj >= 0) {
            data1[jj + 1] = data1[jj];
            jj = jj - 1;
        }//end
        data1[0] = X[ptr];

        ptr = ptr + 1;
    }// end

    int Tpa = (mpa - 1);
    int T = Tpa + 1;
    int j = 0;
    while (j < (dataLength - T)) {
        Xpc[j] = Xpc[j + T];
        j = j + 1;
    }//end

    j = dataLength - (T - 1);
    while (j < dataLength) {
        Xpc[j] = 0;
        j = j + 1;
    }//end


}



//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::GetMinimmumAngle
/// \param leftmost
/// \param rightmost
/// \param

int GetMinimmumAngleB(const float *beat, int leftmost, int rightmost) {

    //  %% look for the minimun value of the angle between two segment
    //  %% having a common mid point and equal lenghts of 10ms
    double minimumAngle = 180.0;
    int temp;
    int QRSstart;

    if (leftmost > rightmost) {
        temp = leftmost;
        leftmost = rightmost;
        rightmost = temp;
    }//end

    if ((leftmost + 2 * BEAT_MS10) > rightmost) {
        // % leftmost=leftmost+BEAT_MS20-1;
        QRSstart = (int) ((leftmost + rightmost) / 2);
        return QRSstart;
    }

    // end
    int leftX;
    // float leftY;
    double leftY;
    int centerX;
    int rightX;
    // float rightY;
    double rightY;
    // float centerY;
    double centerY;
    int temPosition = 0;
    double tempAngle;

    int vectorCenterLeftX;
    // float vectorCenterLeftY;
    double vectorCenterLeftY;
    int vectorCenterRightX;
    double vectorCenterRightY;

    double vectorProductX;
    double vectorA;
    double vectorB;
    double angleValue;


    while ((leftmost + 2 * BEAT_MS10) <= rightmost) {
        leftX = leftmost;
        if (leftX <= 0) {
            // int  kk=0;
            // int  ll=kk;
            leftX = 0;
        }//end
        leftY = beat[leftX];
        centerX = leftmost + BEAT_MS10;
        centerY = beat[centerX];
        rightX = centerX + BEAT_MS10;
        rightY = beat[rightX];

        // %% compute vetor value
        vectorCenterLeftX = leftX - centerX;
        vectorCenterLeftY = leftY - centerY;
        vectorCenterRightX = rightX - centerX;
        vectorCenterRightY = rightY - centerY;


        // %% calculate the vector product
        vectorProductX = vectorCenterLeftX * vectorCenterRightX + vectorCenterLeftY * vectorCenterRightY;
        vectorA = sqrt(vectorCenterLeftX * vectorCenterLeftX + vectorCenterLeftY * vectorCenterLeftY);
        vectorB = sqrt(vectorCenterRightX * vectorCenterRightX + vectorCenterRightY * vectorCenterRightY);
        //  float minimumAngle; minimumAngle   180.0


        // %% looking for the minimum  angle
        angleValue = vectorProductX / (vectorA * vectorB);
        tempAngle = acos(angleValue) * 180 / M_PI;
        if (tempAngle <= minimumAngle) {
            minimumAngle = tempAngle;
            temPosition = centerX;

        }//end

        // %% move to next segment
        leftmost = centerX;
        centerX = rightX;
        rightX = rightX + BEAT_MS10;

    }// end  while((leftmost+2*BEAT_MS10)<=rightmost)

    if (tempAngle >= 179.9999)
        QRSstart = rightmost;
    // return;
    // end
    // %% return the result
    if (temPosition <= rightmost)
        QRSstart = temPosition;
    // end


    return QRSstart;


}








//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::IsoCheckDiffRight
/// \param data
/// \param j
/// \param isoLength


int IsoCheckDiffRightB(float *data, int j, int isoLength, float limit) {

    // %   IsoCheck determines whether the amplitudes of a run
    // %   of data fall within a sufficiently sall amplitude that
    // %	the run can be considered isoelectric.

    int i = j;
    int k = 0;
    // int returnValue=0;
    float max = data[i];
    float min = data[i];
    while (k < isoLength) {
        //   %  if(abs(data(i-k)-data(i-k+1))>=limit)
        if (fabsf(data[i + k] - data[i + k - 1]) >= limit / 2)
            // returnValue=1;
            return 1;
        // end

        if (data[i + k] > max)
            max = data[i + k];
        else if (data[i + k] < min)
            min = data[i + k];
        //end
        // end

        k = k + 1;
    }// end

    if (((fabsf(data[i] - data[i + isoLength - 1])) > 3 * limit))//%%||abs(max-min)>=limit)
    {
        return 1;
    } else {
        return 0;
    }

}


//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::IsoCheckDiff
/// \param data
/// \param j
/// \param isoLength

int IsoCheckDiffB(float *data, int j, int isoLength, float limit) {

    // %   IsoCheck determines whether the amplitudes of a run
    // %   of data fall within a sufficiently small amplitude that
    // %	the run can be considered isoelectric.
    //   % BEAT_MS20=10;


    int i = j;
    int k = 0;
    // %%% search in the interval  from the biggest peak to 120ms
    while ((k < isoLength) && (fabsf(data[i + k] - data[i + k - 1]) < limit))
        k = k + 1;


    if ((k >= isoLength) && ((fabsf(data[i] - data[i + isoLength - 1])) < 3 * limit)) {
        return 0;
    } else {
        return 1;
    }


}


//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::IsoCheck
/// \param data
/// \param j
/// \param isoLength
int IsoCheckB(const float *data, int j, int isoLength) {

    //  BEAT_MS1000=500;
    //  ISO_LIMT=20;
    float max = data[j];
    float min = data[j];
    int i = j;
    int k = 0;
    int returnValue = 0;
    while (k < (isoLength) && ((i + k) < BEAT_MS1000)) {
        if (data[i + k] > max)
            max = data[i + k];
        else if (data[i + k] < min)
            min = data[i + k];
        //end
        // end
        k = k + 1;
    }// end

    if ((max - min) < ISO_LIMT)
        returnValue = 1;
    // return  ;
    // end
    returnValue = 0;


    return returnValue;

}


//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::qrsBoundary
/// \param beat
/// \param Crit
/// \param mainWaveHeight
/// \param QRSbegin
/// \param QRSfinish
/// \param QRSPeaksPositions
/// \param QRSPeaks
/// \param qrsIndex
///
///
///
///
void qrsBoundaryB(const float *beat, float Crit, float mainWaveHeight, int *QRSbegin, int *QRSfinish,
                  int *QRSPeaksPositions,
             float *QRSPeaks, int *qrsIndex) {
    float firstDiff[arrayLimit];
    memset(firstDiff, 0, sizeof(float) * arrayLimit);
    float secondDiff[arrayLimit];
    memset(secondDiff, 0, sizeof(float) * arrayLimit);
    float fnew[arrayLimit];
    memset(fnew, 0, sizeof(float) * arrayLimit);
    float data[arrayLimit];
    memset(data, 0, sizeof(float) * arrayLimit);
    float result[arrayLimit];
    memset(result, 0, sizeof(float) * arrayLimit);
    float y[arrayLimit + 5];
    memset(y, 0, sizeof(float) * (arrayLimit + 5));

    // static int noiseFlagArray[3];
    // static int noiseTag;

    float sum = 0;
    int j = 2;

    int kk = 0;
    while (kk < (arrayLimit + 5)) {
        y[kk] = beat[kk];
        kk++;
    }

    j = 1;
    while (j <= (arrayLimit - 1)) {
        firstDiff[j] = (beat[j + 1] - beat[j - 1]) / 2;
        j = j + 1;
    }//end

    //     fnew=zeros(1,arrayLimit);
    j = 2;
    while ((j > 1) && (j < (arrayLimit - 2))) {
        secondDiff[j] = ((-1) * beat[j - 2] - 2 * beat[j - 1] + 2 * beat[j + 1] + beat[j + 2]) / 8;
        j = j + 1;
    }//end
    j = 1;
    while (j <= arrayLimit) {
        fnew[j] = firstDiff[j] * firstDiff[j] + secondDiff[j] * secondDiff[j];
        j = j + 1;
    }//end

    // data=zeros(1,arrayLimit);
    // result=zeros(1,arrayLimit);

    int ptr = 1;
    //    float peakPosition=0;
    float qrsMaxValue = 0;
    j = 0;
    while (j < arrayLimit) {
        sum = sum + fnew[j];
        sum = sum - data[ptr];
        data[ptr] = fnew[j];
        ptr = ptr + 1;
        if (ptr == BEAT_MS250 + 1)
            ptr = 1;
        //end
        result[j] = sum;

        if (((sum / BEAT_MS250) + 1) > qrsMaxValue)
            qrsMaxValue = sum / BEAT_MS250;
        //end

        j = j + 1;

    }//end  while(j<=arrayLimit)

    //  y=result/BEAT_MS250;

    j = 0;
    while (j < arrayLimit) {
        result[j] = result[j] / BEAT_MS250;
        j++;
    }


    int jj = FIDMARK;
    while ((jj > 2) && (fabsf(result[jj] - result[jj - 2]) < Crit / 2)) //%% for flat top
    {
        jj = jj - 1;
        if (fabsf(result[jj]) < Crit / 10)
            break;
        //end
    }//end
    while ((jj > 6) && ((result[jj] - result[jj - 6]) > (BEAT_MS250 / 60.0)) && (result[jj] > (BEAT_MS250 / 60.0))) {
        jj = jj - 1;
    }//end



    *QRSbegin = jj;
    //  QRSPeaks=zeros(1,BEAT_MS20);
    //  QRSPeaksPositions=zeros(1,BEAT_MS20);

//         kk=0;
//        while(kk<BEAT_MS20)
//        {
//           QRSPeaks[kk]=0;
//           QRSPeaksPositions[kk]=0;
//          kk++;
//        }

    int index = 0;
    *qrsIndex = 0;
    jj = FIDMARK;
    while (((result[jj]) < (qrsMaxValue - 1)) && (jj < arrayLimit))
        jj = jj + 1;
    //end
    *QRSfinish = jj;
    float temp;
    //    y=beat;

    int i = (FIDMARK - BEAT_MS90 - 1);
    int peakIndex = 0;
    float lastTemp = 0;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while ((index < 8) && (i < (BEAT_MS500 + 15)) && ((i + BEAT_MS10) < arrayLimit)) {


        if ((i < (FIDMARK + 5)) && (i > 8)) {
            if ((fabsf(y[i] - y[i - BEAT_MS10]) > Crit) && (fabsf(y[i] - y[i + BEAT_MS10]) > Crit)) {
                if (((((y[i] - y[i - BEAT_MS10])) > 0) && ((y[i] - y[i + BEAT_MS10]) > 0)) ||
                    (((y[i] - y[i - BEAT_MS10]) < 0) && (y[i] - y[i + BEAT_MS10] < 0))) {
                    index = index + 1;
                    kk = i;
                    //  downward
                    if (y[kk] < (y[kk - BEAT_MS10] + y[kk + BEAT_MS10]) / 2) {
                        jj = kk - BEAT_MS10;
                        temp = y[jj];
                        peakIndex = jj;
                        jj = jj + 1;
                        while (jj < (kk + BEAT_MS10)) {
                            if (temp > y[jj]) {
                                temp = y[jj];
                                peakIndex = jj;
                            }
                            // end
                            jj = jj + 1;
                        }//end
                    } else                          //   upward
                    {
                        if ((y[kk] > (y[kk - BEAT_MS10] + y[kk + BEAT_MS10]) / 2) && (y[kk] > y[kk - 2 * BEAT_MS10])) {
                            jj = kk - BEAT_MS10;
                            temp = y[jj];
                            peakIndex = jj;
                            jj = jj + 1;
                            while (jj < (kk + BEAT_MS10)) {
                                if (temp < y[jj]) {
                                    temp = y[jj];
                                    peakIndex = jj;
                                }//end
                                jj = jj + 1;
                            }// end
                        } // if((y[kk]>(y[kk-BEAT_MS10]+y[kk+BEAT_MS10])/2)&&(y[kk]>y[kk-2*BEAT_MS10]))
                    }// end    else
                    // %%%%%%%%%%%%%%%%%% mainWaveHeight
                    if ((peakIndex >= 1) && (fabsf(y[peakIndex]) > mainWaveHeight * 1 / 10))//% 22:38 28/07 2018
                    {
                        if (lastTemp != temp) {

                            QRSPeaksPositions[index - 1] = peakIndex;
                            QRSPeaks[index - 1] = temp;
                            lastTemp = temp;
                            //index=index+1;

                        } else {
                            if (index >= 1) {
                                index--;
                            }
                        }
                    } else {
                        index = index - 1;
                    }
                    i = i + 2 * BEAT_MS6;

                }//end   if((abs(y[i]-y[i-BEAT_MS10])>C
            } // end if((i<(FIDMARK+5))&&(i>8))
        } //end  if((i<(FIDMARK+5))&&(i>8))  592
        else //%%%%%%%%%
            //  %%%% i>FIDMARK
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        {
            if ((i > 4) && (i >= (FIDMARK - BEAT_MS40)) && (i <= (FIDMARK + BEAT_MS40))) {
                if ((fabsf(y[i] - y[i - BEAT_MS8]) > Crit / 2) && (fabsf(y[i] - y[i + BEAT_MS8]) > Crit / 2)) {
                    if (((((y[i] - y[i - BEAT_MS8])) > 0) && ((y[i] - y[i + BEAT_MS8]) > 0)) ||
                        (((y[i] - y[i - BEAT_MS8]) < 0) && (y[i] - y[i + BEAT_MS8] < 0))) {
                        index = index + 1;
                        // downward
                        kk = i;
                        if (y[kk] < (y[kk - BEAT_MS8] + y[kk + BEAT_MS8]) / 2) {
                            jj = kk - BEAT_MS8;
                            temp = y[jj];
                            peakIndex = jj;
                            jj = jj + 1;
                            while (jj < (kk + BEAT_MS8)) {
                                if (temp > y[jj]) {
                                    temp = y[jj];
                                    peakIndex = jj;
                                }
                                //end
                                jj = jj + 1;
                            }// end
                        } else
                            //   upward
                        {
                            if ((y[kk] > (y[kk - BEAT_MS8] + y[kk + BEAT_MS8]) / 2) &&
                                (y[kk] > y[kk + (int) (1.5 * BEAT_MS20)])) // %% avoid little ripple
                            {
                                jj = kk - BEAT_MS6;
                                temp = y[jj];
                                peakIndex = jj;
                                jj = jj + 1;
                                while (jj < (kk + BEAT_MS8)) {
                                    if (temp < y[jj]) {
                                        temp = y[jj];
                                        peakIndex = jj;
                                    }
                                    //end
                                    jj = jj + 1;
                                }
                                //end
                            }//end

                        }//  end
                        //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        // int ll=peakIndex;
                        if ((index > 1) && (fabsf(y[peakIndex]) >= mainWaveHeight * 1 / 3) &&
                            (peakIndex != (QRSPeaksPositions[index - 1])))//%0.58)
                        {
                            if (lastTemp != temp) {
                                QRSPeaksPositions[index - 1] = peakIndex;
                                QRSPeaks[index - 1] = temp;
                                lastTemp = temp;

                            } else {
                                if (index >= 1) {
                                    index--;
                                }
                            }
                        } else {
                            index = index - 1;
                        }
                        //end
                        //  %
                        i = i + 2 * BEAT_MS10;
                    } //end
                } // end
            } else      //%%%%%%%% end for  if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
                //    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            {
                if ((i > BEAT_MS10) && (fabsf(y[i] - y[i - BEAT_MS10]) > Crit / 2) &&
                    (fabsf(y[i] - y[i + BEAT_MS10]) > Crit / 2)) {
                    if ((((y[i] - y[i - BEAT_MS10]) > 0) && ((y[i] - y[i + BEAT_MS10])) > 0) ||
                        (((y[i] - y[i - BEAT_MS10]) < 0) && (y[i] - y[i + BEAT_MS10] < 0))) {
                        index = index + 1;
                        kk = i;

                        //  downward
                        if ((y[kk] < (y[kk - 2 * BEAT_MS10] + y[kk + 2 * BEAT_MS10]) / 2) &&
                            (y[kk] < y[kk - 2 * BEAT_MS20])) {
                            jj = kk - BEAT_MS10;
                            temp = y[jj];
                            peakIndex = jj;
                            jj = jj + 1;
                            while (jj < (kk + BEAT_MS10)) {
                                if (temp > y[jj]) {
                                    temp = y[jj];
                                    peakIndex = jj;
                                } //end
                                jj = jj + 1;
                            } // end
                        } else
                            //   upward
                        {
                            if ((y[kk] > (y[kk - 2 * BEAT_MS10] + y[kk + 2 * BEAT_MS10]) / 2) &&
                                (y[kk] > y[kk + 2 * BEAT_MS20]))//  %% avoid little ripple
                            {
                                jj = kk - BEAT_MS6;
                                temp = y[jj];
                                peakIndex = jj;
                                jj = jj + 1;
                                while (jj < (kk + BEAT_MS10)) {
                                    if (temp < y[jj]) {
                                        temp = y[jj];
                                        peakIndex = jj;
                                    }//end
                                    jj = jj + 1;
                                }//end
                            }//end
                        }//end

                        //   %  if(abs(y(peakIndex))>abs(y(FIDMARK)/2))
                        if ((index > 1) && (fabsf(y[peakIndex]) >= mainWaveHeight * 1 / 4) &&
                            (peakIndex != QRSPeaksPositions[index - 1]))//%0.58)  22:38  28/07 2018
                        {
                            if (lastTemp != temp) {

                                QRSPeaksPositions[index - 1] = peakIndex;
                                QRSPeaks[index - 1] = temp;
                                lastTemp = temp;
                                // index=index+1;

                            } else {
                                if (index >= 1) {
                                    index--;
                                }
                            }
                        } else {
                            index = index - 1;
                        }

                        i = i + 2 * BEAT_MS10;
                    }// end  729
                }//end  line 726    if((i>BEAT_MS10)&&(fabs(y[i]-y[i-BEAT_MS10])>Crit/2)&&(fabs(y[i]-y[i+BEAT_MS10])>Crit/2))

            }//end  line 724

        } //end    line 563  if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))

        i = i + 1;

        *qrsIndex = index;


    } // end  while ((i<(BEAT_MS500+15))&&((i+BEAT_MS10)<arrayLimit))


}



//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::BeatFeaturesGet
/// \param tempBeat
/// \param Crit
/// \param mainWaveHeight
/// \param QRSbegin
/// \param QRSfinish
/// \param QRSPeaksPositions
/// \param QRSPeaks
/// \param index
///
///
///
///

void BeatFeaturesGet(BeatFeaturesGet_CTX *ctx, float *beat, int *SPeak, int interval, int *onset,
                     int *offset, int *QPeak,               float *speakheight,                float *rweight, int *noiseFlag, float *isoelectricLevel, int *notchFlag, int init) {

//        *isoelectricLevel=0;
//        *QPeak=180;
//        *onset=180;

    float QRSPeaks[10];//=zeros(1,BEAT_MS20);
    memset(QRSPeaks, 0, sizeof(float) * 10);
    int QRSPeaksPositions[10];//=zeros(1,BEAT_MS20);
    memset(QRSPeaksPositions, 0, sizeof(float) * 10);
    float Xpc[sampleRate];//=zeros(1,dataLength);


    float *noiseFlagArray = ctx->noiseFlagArray;
    int *ctx_noiseTag = &(ctx->noiseTag);


    if (init) {
        //noiseFlagArray[3]={0};
        memset(noiseFlagArray, 0, sizeof(float) * 3);
        (*ctx_noiseTag) = 0;
        return;
    }

    int kk;

    int i = (int) (FIDMARK - BEAT_MS80);
    float maxQRS = beat[FIDMARK];
    float minQRS = beat[FIDMARK];
    int maxPosition = i;
    int minPosition = i;
    while (i < (FIDMARK + BEAT_MS160)) {
        if ((beat[i] > maxQRS) && ((beat[i] > beat[i + 4]) && (beat[i] > beat[i - 4]))) {
            maxQRS = beat[i];
            maxPosition = i;
        } else if (beat[i] < minQRS) {
            minQRS = beat[i];
            minPosition = i;
        }//end
        // end

        i = i + 1;
    }//end
    float Crit = 0.02f * (maxQRS - minQRS);


//        %find S
    int sWavePosition = 0;
    int sWaveFlag = 0;
    int start = (minPosition - BEAT_MS10); //%default
    if (start < FIDMARK)
        start = FIDMARK + 1;
    // end
    int nextqrsindex = (minPosition);
    int sgn = (int) (beat[start] - beat[nextqrsindex]);
    if (sgn != 0)
        sgn = sgn / abs(sgn);
    else
        sgn = 1;
    //end
    float currentSignal;
    float nextSignal;
    float newsign;
    int jj = start + 1;
    while (jj < (FIDMARK + BEAT_MS160)) {
        currentSignal = beat[jj];
        nextSignal = beat[jj + 1];
        newsign = currentSignal - nextSignal;
        if (newsign != 0)
            newsign = newsign / fabsf(newsign);
        else
            newsign = 1;
        // end
        if (newsign != (float) sgn) {
            sWavePosition = jj;
            jj = (FIDMARK + BEAT_MS160);// %safe gaurd for the break statement
            break;
        }
        // end
        jj = jj + 1;
    }//end
    int speakTrueFlag = 0;
//          if(sWavePosition<(FIDMARK+ BEAT_MS160))
//             sWaveFlag=1;
//      // %       SPeak=sWavePosition;
//      // %       speakTrueFlag=1;
//           x=0;
//          end
    int mFrom = minPosition - 2 * BEAT_MS10;
    if (mFrom < 0)
        mFrom = 0;
    // end

    int mTo = minPosition + 2 * BEAT_MS10;
    if (mTo > (FIDMARK + BEAT_MS160))
        mTo = FIDMARK + BEAT_MS160;
    //end
    //%       minSlopeLeft=beat(minPosition)-beat(mFrom);
    // %       minSlopeRight=beat(minPosition)-beat(mTo);
    float middleValue = (beat[mFrom] + beat[mTo]) / 2;
    if (beat[minPosition] < middleValue)
        sWaveFlag = 1;
    // end


    float mainWaveHeight;
    if ((fabsf(maxQRS)) > fabsf(minQRS))
        mainWaveHeight = fabsf(maxQRS);
    else
        mainWaveHeight = fabsf(minQRS);
    // end
/////////////////////////////////////////////////////////////////
    //  utility tempPointer;
    float mnoise0 = HFnoiseCheck(beat);
    float mnoise1 = noiseLevel(beat, 500);
    if (mnoise0 > 60 || mnoise1 > 60) {
        mnoise0 = 0;
        *noiseFlag = 1;
        return;
    }

//////////////////////////////////////////////////////////////////
    // [QRSbegin ,QRSfinish,QRSPeaksPositions,QRSPeaks,qrsIndex]= QRSBoundary( beat ,Crit,mainWaveHeight);
    //int*QRSbegin=new int[1];
    // int *QRSfinish=new int[1];
//              int  *qrsIndex=new int[1];
    int QRSbegin = 0;
    int QRSfinish = 0;
    int qrsIndex = 0;

    //   void qrsBoundary(float*beat ,float Crit,float mainWaveHeight,int*QRSbegin ,int*QRSfinish,int*QRSPeaksPositions,float*QRSPeaks,int* qrsIndex);
    qrsBoundaryB(beat, Crit, mainWaveHeight, &QRSbegin, &QRSfinish, QRSPeaksPositions, QRSPeaks, &qrsIndex);

    if (qrsIndex <= 0) {
        qrsIndex = 0;
        QRSPeaksPositions[qrsIndex] = FIDMARK;
        if (beat[FIDMARK] < 0)
            QRSPeaks[qrsIndex] = minQRS;
        else
            QRSPeaks[qrsIndex] = maxQRS;
        // end
    }
    // noiseFlag=0;

    // dataLength=length(beat);
    //  sampleRate=500;
    // x=0:1/sampleRate:dataLength/sampleRate-1/sampleRate;
    //%25;%50;%100;%50;
    // X=beat;
    //[ returnValue ] = HFNoiseCheck(X);
    // X=beat;
    // X  = LowPassFilter(X,Fs,Fpa);
    // Xaux=X(80:170);
    // mnoise=noiselevel(Xaux);

    int Fpa = 50;
    int start1 = 80;
    int finish = 170;
    int len = finish - start1;
    float XX[500] = {0};


    float Xaux[90];
    float mnoise = 0;
    //      utility tempPointer;
    float HFlevel = HFnoiseCheck(beat);
    LowPassFilterB(beat, Fs, Fpa, XX);
    ecgCopy(XX, Xaux, start1, finish);
    //mnoise=tempPointer.noiseLevel(Xaux,len);
    mnoise = noiseLevel(Xaux, len);
    int leng = 90;
    // float xMax=tempPointer.ECGMax( Xaux,leng );
    // float xMin= tempPointer.ECGMin( Xaux,leng );
    // float motionSum=tempPointer.ecgSum(Xaux,leng);
    // pp=motionSum/leng;
    start1 = 70;
    finish = 160;
    len = finish - start1;
    ecgCopy(XX, Xaux, start1, finish);
    // Xaux=X(70:160);
    float artifactsIndex = ecgSum(Xaux, leng);


    kk = 0;
    while (kk < 2) {
        noiseFlagArray[kk] = noiseFlagArray[kk + 1];
        kk = kk + 1;
    }
    // end
    noiseFlagArray[2] = (*ctx_noiseTag);
    if (((noiseFlagArray[0] == 1) && (noiseFlagArray[1] == 1)) ||
        ((noiseFlagArray[1] == 1) && (noiseFlagArray[2] == 1)) ||
        ((noiseFlagArray[0] == 1) && (noiseFlagArray[2] == 1)))
        if (mnoise > 9) {
            *noiseFlag = 1;
            mnoise = 33;
        }//end
    // end

    if ((mnoise > 18) && (mnoise < 32) && (fabsf(minQRS) > maxQRS))
        mnoise = 0;
    //end

    //   %  if((returnValue >hfThreshold)||(mnoise>35))%% 45
    if ((HFlevel > 75) || (mnoise > 24) ||
        ((artifactsIndex > 80) && (mnoise > 19.1)))//    % Change from 61  46 to 51  35
    {
        *noiseFlag = 1;
        if (mnoise > 18) {
            (*ctx_noiseTag) = 1;
            return;
        }
    } else {
        (*ctx_noiseTag) = 0;
    }

    int qrsType = 0;
    int sPeakFlag = 0;
    int qPeakFlag = 0;
//                   //QPeak=0;
//                   //SPeak=0;
    int pqNoise = 0;
//                  // notchFlag=0;
    sWaveFlag = 0;

    int muitiPeakFlag = 0;
    // if(qrsIndex>=1)
    // {
    switch (qrsIndex) {

        case 1:
            if ((QRSPeaks[0] > 0) && (minQRS < 0) && (minPosition > FIDMARK) && (minPosition <= 260)) {
                sPeakFlag = TRUE;
                *SPeak = minPosition;
                sWaveFlag = TRUE;
            } else if ((qrsIndex == 1) && (QRSPeaks[0] < 0)) {
                sPeakFlag = TRUE;
                *SPeak = QRSPeaksPositions[0];//minPosition;
                sWaveFlag = TRUE;
                speakTrueFlag = TRUE;
            } // end
            //end
            QRSPeaksPositions[0] = FIDMARK;
            break;
        case 2:
            if ((QRSPeaksPositions[0] <= (FIDMARK + 2) && (QRSPeaksPositions[0] >= (FIDMARK - 2))) &&
                (QRSPeaks[1] < 0)) {
                sPeakFlag = TRUE;
                *SPeak = QRSPeaksPositions[1];
                speakTrueFlag = TRUE;
            } else if (QRSPeaksPositions[1] <= (FIDMARK + 3)) {
                if (QRSPeaks[0] < 0) {
                    qPeakFlag = TRUE;
                    *QPeak = QRSPeaksPositions[0];
                } //end
                if ((QRSPeaks[0] > 0) && (fabsf(QRSPeaks[1]) > QRSPeaks[0] * 2) &&
                    (QRSPeaksPositions[1] > (FIDMARK + 15))) {
                    sPeakFlag = TRUE;

                    *SPeak = QRSPeaksPositions[1];
                    speakTrueFlag = TRUE;
                } else {
                    if (QRSPeaks[1] < 0) {
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[1];
                        speakTrueFlag = TRUE;
                    }

                }// end   else
            } // end    if( QRSPeaksPositions[1]<=(FIDMARK+3))
            else if ((QRSPeaks[0] < 0) && (QRSPeaks[1] > 0) && (fabsf(QRSPeaks[0]) > fabsf(QRSPeaks[1]))) {
                sPeakFlag = TRUE;
                *SPeak = QRSPeaksPositions[0];//minPosition;
                sWaveFlag = TRUE;
                speakTrueFlag = TRUE;
            }

            // end
            //end

            if ((QRSPeaks[0] > 0) && (QRSPeaks[1] > 0))
                *notchFlag = TRUE;
            // end

            break;

        case 3: {

            if ((QRSPeaks[1] > QRSPeaks[0]) && (QRSPeaks[1] > QRSPeaks[2])) {
                if (beat[QRSPeaksPositions[1]] >
                    ((beat[QRSPeaksPositions[1] - BEAT_MS10] + beat[QRSPeaksPositions[1] + BEAT_MS10])) / 2) {
                    if (beat[QRSPeaksPositions[0]] <
                        ((beat[QRSPeaksPositions[0] - BEAT_MS10] + beat[QRSPeaksPositions[0] + BEAT_MS10])) / 2) {
                        qPeakFlag = TRUE;
                        *QPeak = QRSPeaksPositions[0];
                    }//end
                    if (beat[QRSPeaksPositions[2]] <
                        ((beat[QRSPeaksPositions[2] - BEAT_MS10] + beat[QRSPeaksPositions[2] + BEAT_MS10])) / 2) {
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[2];
                        speakTrueFlag = TRUE;
                    } // end
                }//end
                else {
                    if (QRSPeaks[1] >= 0.95 * mainWaveHeight) {
                        qPeakFlag = TRUE;
                        *QPeak = QRSPeaksPositions[1];
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[2];
                    } else {
                        if (QRSPeaks[2] >= 0.95 * mainWaveHeight) {
                            // %                      qPeakFlag=1;
                            // %                      QPeak=QRSPeaksPositions[1];
                        } else {
                            if ((QRSPeaks[0] >= 0.6 * mainWaveHeight) && (QRSPeaks[1] >= 0.5 * mainWaveHeight) &&
                                (QRSPeaks[2] >= 0.5 * mainWaveHeight)) {
                                qrsType = 3;
                            } else {
                                if (QRSPeaks[0] >= 0.95 * mainWaveHeight) {
                                    sPeakFlag = TRUE;
                                    *SPeak = QRSPeaksPositions[2];
                                    speakTrueFlag = TRUE;
                                }//end
                            } //end
                        }//end

                    }//end
                }// end
            }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRUE

            jj = 0;
            while (jj < qrsIndex) {
                if ((QRSPeaksPositions[jj] > (FIDMARK - BEAT_MS10)) && (QRSPeaksPositions[jj] < (FIDMARK + BEAT_MS10)))
                    break;
                //end
                jj = jj + 1;
            }//end

            switch (jj) {
                case 0:
                    if (QRSPeaks[2] < 0) {
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[2];
                        speakTrueFlag = TRUE;
                    }// end
                    break;
                case 1:
                    if (QRSPeaks[0] < 0) {
                        qPeakFlag = TRUE;
                        *QPeak = QRSPeaksPositions[0];
                    }// end

                    if (QRSPeaks[2] < 0) {
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[2];
                        speakTrueFlag = TRUE;
                    }//end
                    break;
                case 2:
                    if ((QRSPeaks[0] > 0) && (QRSPeaks[1] > 0) && (QRSPeaks[2] > 0))
                        *notchFlag = TRUE;
                    //end
                    if (QRSPeaks[0] < 0) {
                        qPeakFlag = TRUE;
                        *QPeak = QRSPeaksPositions[0];
                    }
                    break;
                default:;
                    //end
            }//end
            //     %%%%% Notch Identification
            if ((maxQRS) > fabsf(minQRS)) {
                if ((QRSPeaks[0] > 0) && QRSPeaks[1] > 0)
                    *notchFlag = TRUE;
                //end
            } else {
                if ((QRSPeaks[0] < 0) && QRSPeaks[1] < 0)
                    *notchFlag = TRUE;
            } //end
            break;
        }  //end   case 3:

// %%%%%%%%%%%%%%%% end case 3 %%%%%%%%%%%


        case 4:
            //    %%%%%% look for FIDMARK %%%%%%%
            //  int peakIndex=0;
            jj = 0;
            while (jj < qrsIndex) {
                if ((QRSPeaksPositions[jj] > FIDMARK - BEAT_MS10) && (QRSPeaksPositions[jj] < FIDMARK + BEAT_MS10))
                    break;
                //end
                jj = jj + 1;
            }//end
            //  %%%%%% Q position  %%%%%%
            if ((jj >= 2) && (qrsIndex >= 2) && (QRSPeaks[0] < 0)) {
                qPeakFlag = TRUE;
                *QPeak = QRSPeaksPositions[1];
            }
            //               else
            //                  if((peakIndex>=2)&&((QRSPeaks[2]<0))&&(fabs(QRSPeaks[2])>QRSPeaks(1)/2))
            //  %                    qPeakFlag=1;
            //  %                    QPeak=QRSPeaksPositions(2);TRUE;
            //                  end
            //                end

            if ((jj == 3) && (QRSPeaks[0] > 0) && (QRSPeaks[1] > 0) && (QRSPeaks[2] > 0) && (QRSPeaks[3] < 0)) {
                *notchFlag = TRUE; // %% segmenttation
                sPeakFlag = TRUE;
                *SPeak = QRSPeaksPositions[3];
                speakTrueFlag = TRUE;
            }//end
            if (((QRSPeaks[0] > 0) && (QRSPeaks[2] > 0)) || ((QRSPeaks[1] > 0) && (QRSPeaks[2] > 0)) ||
                ((QRSPeaks[1] > 0) && (QRSPeaks[3] > 0)))
                *notchFlag = TRUE;  //%% segmenttation

            //end

            break;


        default:

            //     int muitiPeakFlag=0;
            if (beat[FIDMARK] > 0) {
                kk = 0;
                int qrsLogic;
                while ((kk < qrsIndex)) {
                    qrsLogic = (QRSPeaksPositions[kk] > 175) && (QRSPeaksPositions[kk] < BEAT_MS500);
                    if (qrsLogic && (QRSPeaks[kk] < 0) && (fabsf(QRSPeaks[kk]) > mainWaveHeight / 2))
                        muitiPeakFlag = muitiPeakFlag + 1;
                    //end
                    kk = kk + 1;
                }//end
            }


            if (((QRSPeaks[0] > 0) && (QRSPeaks[2] > 0)) || ((QRSPeaks[1] > 0) && (QRSPeaks[2] > 0)) ||
                ((QRSPeaks[1] > 0) && (QRSPeaks[3] > 0)))
                *notchFlag = TRUE; // %% segmenttation
            // end


            if (muitiPeakFlag >= 2) {
                qPeakFlag = TRUE;
                *QPeak = QRSPeaksPositions[0];
                sPeakFlag = TRUE;
                *SPeak = QRSPeaksPositions[qrsIndex];
                speakTrueFlag = TRUE;
            } else
// %%%%%%%%%%% upward %%%%%%%%%%%
            {
                if ((maxQRS) > fabsf(minQRS)) {
                    jj = 0;
                    while ((QRSPeaks[jj] < mainWaveHeight / 2) && (jj < qrsIndex))
                        jj = jj + 1;
                    //end
                    int j = jj - 1;
                    if (j <= 1) {
                        qPeakFlag = 1;
                        *QPeak = QRSPeaksPositions[0];
                    } else if (((QRSPeaksPositions[j]) < (FIDMARK - BEAT_MS10)) && (QRSPeaks[j] < 0)) {
                        pqNoise = j;
                        qPeakFlag = 1;
                        *QPeak = QRSPeaksPositions[j];
                    }
                    //end
                    // end  %%   if(j<=1)
                } else
                    // %%%%%%%%%%% downward %%%%%%%
                {
                    jj = 0;
                    while ((fabsf(QRSPeaks[jj]) < mainWaveHeight / 3) && (jj < qrsIndex))
                        jj = jj + 1;
                    //end
                    int j = jj - 1;
                    //               if(j<=1)
                    // %                  qPeakFlag=1;
                    // %                  QPeak=QRSPeaksPositions(1);
                    //               else
                    // % if((QRSPeaksPositions(j))<(FIDMARK-BEAT_MS10))&&(QRSPeaks(j)<0)
                    if ((j >= 3) && (QRSPeaks[j] < 0))
                        pqNoise = j;
                    //  %                   qPeakFlag=1;
                    //  %                   QPeak=QRSPeaksPositions(1);
                    //                  end
                    //                end   %%% if(j<=1)
                    qPeakFlag = 1;
                    *QPeak = QRSPeaksPositions[1];

                }//  end  %% end    if ((maxQRS )>abs(minQRS))

                jj = 0;
                while ((jj < qrsIndex) && (QRSPeaksPositions[jj] < (FIDMARK + BEAT_MS10)))
                    jj = jj + 1;
                //end
                if (QRSPeaks[jj] < 0) {
                    sPeakFlag = TRUE;
                    *SPeak = QRSPeaksPositions[jj];
                    speakTrueFlag = TRUE;
                } else {
                    if (qrsIndex >= 1) {
                        sPeakFlag = TRUE;
                        *SPeak = QRSPeaksPositions[qrsIndex];
                        speakTrueFlag = TRUE;
                    }
                }
            }
            break;
    }
    //	} // end if
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch %%%%%%%%%%%%%%%%%%%%%

    i = QRSPeaksPositions[0] - BEAT_MS20;//%%FIDMARK-BEAT_MS20 ;
    while ((i < BEAT_MS1000) && (i > FIDMARK - BEAT_MS120) && (IsoCheckDiffB(beat, i, ISO_LENGTH, Crit) == 1))
        i = i - 1;
    //end
    //%% flat top
    int flatFlag = 0;
//                    float ppp=(beat[i]/beat[FIDMARK]);
//                    float kkk=ppp>4/5;
//                    kkk=ppp>0.8;
    if ((i == (FIDMARK - BEAT_MS20)) || ((float) (beat[i] / beat[FIDMARK]) > 0.8)) {
        flatFlag = 1;
        i = FIDMARK - BEAT_MS20 * 2;
        while ((i < BEAT_MS1000) && (i > FIDMARK - BEAT_MS120) && (IsoCheckDiffB(beat, i, BEAT_MS20, Crit) == 1)) {
            i = i - 1;
        }
        //end
    }//end   if((i==(FIDMARK-BEAT_MS20))||(((beat[i]/beat[FIDMARK])>4/5)))

    int leftMost = i;
    while ((i < BEAT_MS1000) && (i > FIDMARK - BEAT_MS120) && (IsoCheckDiffB(beat, i, BEAT_MS16, Crit / 5) == 1))
        i = i - 1;
    // end
    int iosStart = i;

// %%%%%%%%%%%%%%%%%%%%%%%%%%
    kk = QRSPeaksPositions[0] - BEAT_MS16;//%FIDMARK-8;
    float temp = beat[kk];
    int Qindex = 0;
    while (kk <= (QRSPeaksPositions[0] + BEAT_MS10)) {
        if (beat[QRSPeaksPositions[0]] > 0) {
            if (((beat[kk + 1])) > temp) {
                temp = (beat[kk + 1]);
                Qindex = kk + 1;
            }//end
        } else {
            if ((fabsf(beat[kk + 1])) > temp) {
                temp = fabsf(beat[kk + 1]);
                Qindex = kk + 1;
            }//end
        }//end  else

        kk = kk + 1;
    }//end

    int RPeak = Qindex;


    int j = RPeak - BEAT_MS10;
    while ((j > leftMost) && ((qPeakFlag == 0))) {
        if ((fabsf(beat[j] - beat[j - BEAT_MS10]) > Crit / 2) && (fabsf(beat[j] - beat[j + BEAT_MS10]) > Crit / 2))
            if (((((beat[j] - beat[j - BEAT_MS10])) > 0) && ((beat[j] - beat[j + BEAT_MS10]) > 0)) ||
                (((beat[j] - beat[j - BEAT_MS10]) < 0) && (beat[j] - beat[j + BEAT_MS10] < 0))) {
                if (beat[j] < 0) {
                    *QPeak = j;
                    j = leftMost - BEAT_MS30;//%15;
                    break;
                }
            }//end

        //end
        j = j - 1;
    }//end

    // %%%%%%%%%%%%%%%%%%%%%%%%%%

    if ((*QPeak != 0) && (qPeakFlag == 0)) //EPSILON
        //    if((fabs(*QPeak)<EPSILON)&&(qPeakFlag==0))
    {
        j = *QPeak - BEAT_MS10;
        temp = (beat[j]);
        Qindex = 0;
        while (j <= (*QPeak + BEAT_MS10)) {
            if (beat[FIDMARK] > 0) {
                if (((beat[j + 1])) < temp) {
                    temp = (beat[j + 1]);
                    Qindex = j + 1;
                }//end
            } else if (((beat[j + 1])) > temp) {
                temp = (beat[j + 1]);
                Qindex = j + 1;
            }//end

            // end
            j = j + 1;
        } //end
        if (*QPeak < 0) {
            *QPeak = Qindex;
        }
    }// end


    //%%%%%%%%%%%%%%%%%%%%
    //%% flat top

    if (flatFlag == 1) {
        i = FIDMARK;
        while ((i > FIDMARK - BEAT_MS120) && (fabsf(beat[i] - beat[i - 1]) < Crit))
            i = i - 1;
        //end
        j = i;
    } else
        j = QRSPeaksPositions[0];
    // end
    //  %%%%%%%% for Q slope %%%%%%%%%%%%%%%%%
    int QSlope = 0;
    int qFlag = 0;
    int k;
    while (j > (QRSPeaksPositions[0] - BEAT_MS120)) {
        k = 0;
        int positive = 0;
        int negative = 0;

        while (k < slopeLength) {
            if (((beat[j - k - 2] - beat[j - k]) > 0) && (beat[j - k - 2] - beat[j - k]) > Crit / 2) {
                positive = positive + 1;
                if (positive == slopeLength - 1) {
                    QSlope = j - slopeLength / 2;
                    j = BEAT_MS120;
                    k = slopeLength;
                    qFlag = 1;
                    break;
                }//end
            } else {
                if (((beat[j - k - 2] - beat[j - k]) < 0) && (beat[j - k - 2] - beat[j - k]) < Crit / 2) {
                    negative = negative + 1;
                    if (negative == slopeLength - 2)   //%%%%%%%%%%%%%%%%%%%%%%% 09:58 17/01 2018
                    {
                        QSlope = j - slopeLength / 2;
                        j = BEAT_MS120;
                        k = slopeLength;
                        qFlag = 1;
                        break;
                    }//end
                }//end
            }//end
            k = k + 1;
        }//end

        if (qFlag == 1)
            break;
        //end
        j = j - 1;
    }//end
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ((*QPeak == 0) && (QSlope == 0))
        *QPeak = QRSbegin; // %%%%%%%%%  14:15 12/02 2018 ypz
    //end
    int QRSstart = 0;
    if (*QPeak != 0)
        QRSstart = GetMinimmumAngleB(beat, leftMost, *QPeak);
    else {
        if (QSlope != 0) {
            if (leftMost > QSlope) {
                QRSstart = leftMost;
            } else {
                QRSstart = GetMinimmumAngleB(beat, leftMost, QSlope);
            }//end
        }// end
    }// end

    if ((pqNoise >= 2) && (qPeakFlag == 1)) {
        if ((QRSbegin > QRSstart) && (QRSbegin - QRSstart) > 10) {
            QRSstart = QRSbegin - 10;
        }//end
    }//end



// %%%%%%%%%%%%%%%%%%%%%% QRSstart end   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // %%%%%%%%%%%%%%%%%%%%%% SPeak start  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ((qrsIndex <= 2) && (qrsIndex >= 1) || qrsIndex <= 0)
        j = FIDMARK + BEAT_MS12;
    else
        j = QRSPeaksPositions[qrsIndex - 1] + BEAT_MS6;


    int min = (int) beat[j];
    int SpeakIndex = j;

    // %   while((sPeakFlag==0)&&(j<(QRSPeaksPositions(qrsIndex)+BEAT_MS160)))
    // %   speakTrueFlag
    while ((speakTrueFlag == 0) && (qrsIndex >= 1) && (j < (QRSPeaksPositions[qrsIndex - 1] + BEAT_MS160))) {
        if (beat[j] < (float) min) {
            min = (int) beat[j];
            SpeakIndex = j;
        }

        if ((fabsf(beat[j] - beat[j - BEAT_MS10]) > Crit) && (fabsf(beat[j] - beat[j + BEAT_MS10]) > Crit)) {
            if (((((beat[j] - beat[j - BEAT_MS10])) > 0) && ((beat[j] - beat[j + BEAT_MS10]) > 0)) ||
                (((beat[j] - beat[j - BEAT_MS10]) < 0) && (beat[j] - beat[j + BEAT_MS10] < 0))) {
                if ((j < BEAT_MS520) && (j > FIDMARK) && (fabsf(beat[j]) > Crit) && (beat[j] < 0)) {
                    //  if((j<235)&&(j>FIDMARK)&&(fabs(beat[j])>=mainWaveHeight*1/4))
                    {
                        *SPeak = j;
                        speakTrueFlag = 1;
                        j = FIDMARK + BEAT_MS160 + 1;// %%%%%%%%  Break from loop routine speakThrehold
                        break;
                    }
                }
            }
        }
        j = j + 1;
    }

    // %%%%%%%%%%%%%%%%%%%%%% SPeak end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //  %%%%%%%%%%%%%%%%%%%%%%%%% 10:50 14/02 2018 ypz

    int peakFlag = 0;
    if ((*SPeak == 0) && (sPeakFlag == 0) && (qrsType == 0))
        if ((SpeakIndex != 0) && (SpeakIndex < 280))//%||(SpeakIndex-2~=0))
            if ((beat[FIDMARK] > 0) && (beat[SpeakIndex] < 0) &&
                ((beat[SpeakIndex] - beat[SpeakIndex - 2]) * (beat[SpeakIndex + 2] - beat[SpeakIndex]) < 0)) {
                *SPeak = SpeakIndex;
                peakFlag = 1;
            }


    if ((*SPeak != 0) && (sPeakFlag == 0)) {
        j = *SPeak - BEAT_MS10;
        temp = beat[j];
        int Sindex = 0;
        while (j <= (*SPeak + BEAT_MS10)) {
            if (beat[FIDMARK] > 0) {
                if (((beat[j + 1])) < temp) {
                    temp = (beat[j + 1]);
                    Sindex = j + 1;
                }
            } else {
                if (((beat[j + 1])) > temp) {
                    temp = beat[j + 1];
                    Sindex = j + 1;
                }

            }
            j = j + 1;
        }
        *SPeak = Sindex;
    }

//  %%%%%%%%%% for S slope %%%%%%%%%%%%%%

    int SSlope = 0;
    int sFlag = 0;
    int positive;
    int negative;

    if (qrsIndex <= 1)
        j = QRSPeaksPositions[0] + BEAT_MS20 + 2;
    else {
        if ((qrsIndex - 1) >= 0) {
            j = QRSPeaksPositions[qrsIndex - 1] + BEAT_MS10;
        }
    }
    // end
    // %   while(j<(QRSPeaksPositions(1)+BEAT_MS160))  FIDMARK
    while (j < (FIDMARK + BEAT_MS160)) {
        k = 0;
        positive = FALSE;
        negative = FALSE;

        while (k <= slopeLength) {

            if (((beat[j - k - 3] - beat[j - k]) > 0) && ((beat[j - k - 3]) - (beat[j - k])) > Crit / 2) {
                positive = positive + 1;
                negative = FALSE;
                if ((positive == slopeLength) && (beat[j - slopeLength / 2] < beat[FIDMARK] * 0.8)) {
                    SSlope = j - slopeLength / 2;
                    j = FIDMARK + BEAT_MS160;
                    sFlag = 1;
                    break;
                }
            } else {
                if (((beat[j - k - 3] - beat[j - k]) < 0) && ((beat[j - k] - (beat[j - k - 3]))) > Crit) {
                    negative = negative + 1;
                    positive = FALSE;
                    if ((negative == slopeLength) &&
                        (beat[j - slopeLength / 2] < beat[FIDMARK] * 0.8))//%% for notched peak
                    {
                        SSlope = j - slopeLength / 2;
                        j = FIDMARK + BEAT_MS160;
                        sFlag = 1;
                        break;
                    }
                }
            }
            k = k + 1;
        }
        if (sFlag == 1)
            break;
        //end
        j = j + 1;
    }
// %%%%%%%%%%%%%%%
    int QRSIso;
    int QRSIsoSecond;
    if ((qrsIndex) >= 1) {
        i = QRSPeaksPositions[qrsIndex - 1] + BEAT_MS30;
    }
    if ((*SPeak != 0) && (*SPeak < speakThrehold))
        i = *SPeak;
    //end
    if ((qrsIndex == 1) && (QRSPeaks[qrsIndex] < 0))
        i = QRSPeaksPositions[qrsIndex] + BEAT_MS30;
    //end
    //%% if((speakTrueFlag==0)&&(SPeak~=0)&&(QRSPeaks(1)>0)&&(minQRS<0)&&(minPosition>FIDMARK)&&(minPosition<=260))%%  &&(abs(beat(SPeak))>mainWaveHeight/5)
    if ((speakTrueFlag == 0) && (*SPeak != 0) && (QRSPeaks[0] > 0) && (minQRS < 0) &&
        (fabsf(beat[(int) (*SPeak)]) > mainWaveHeight / 5) && (minPosition <= 260))  //%%%  14:38 08/01 2019 ypz
        i = (int) ((float) *SPeak + BEAT_MS20 * 1.5f);


    while ((i < BEAT_MS1000) && (i < FIDMARK + BEAT_MS400) && (IsoCheckDiffRight(beat, i, ISO_LENGTH_DIFF, Crit) == 1))
        i = i + 1;

    if (speakTrueFlag == 1)
        QRSIso = (int) (i +
                        ISO_LENGTH_DIFF * 1.5);//%isoEnd;% i+ISO_LENGTH_DIFF/2;%%%%%%%%%%%%%%%%%%%%%  10:04 30/10 2017
    else
        QRSIso = (i + BEAT_MS20);


    if (speakTrueFlag == 0) {
        i = FIDMARK + BEAT_MS20;
        while ((i < (FIDMARK + BEAT_MS400)) && (IsoCheckDiffRight(beat, i, ISO_LENGTH_DIFF, Crit) == 1))
            i = i + 1;
        //end
        QRSIsoSecond = (i + BEAT_MS20);

    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15:38 21/12 2018 Ypz
    int stJPoint = 0;
    if ((qrsIndex - 1) >= 0) {
        int isoIndex = QRSPeaksPositions[qrsIndex - 1] + BEAT_MS30;
        while (((isoIndex < FIDMARK + BEAT_MS400)) && (IsoCheckDiffRight(beat, isoIndex, ISO_LENGTH_DIFF, Crit) == 1))
            isoIndex = isoIndex + 1;

        if (SSlope > 1)
            stJPoint = GetMinimmumAngleB(beat, SSlope, isoIndex + BEAT_MS20);
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int localMax = (int) fabsf(beat[280]);
    int localMaxPosition = 281;
    // int pacedFlag=FIDMARK-QRSstart;
    k = 271;
    if (QRSIso > 280)//%%&&(pacedFlag>=BEAT_MS24)
    {
        while (k < 350) {
            if (fabsf(beat[k]) > (float) localMax) {
                //% if((beat(k))>localMax)
                localMax = (int) fabsf(beat[k]);
                localMaxPosition = k;
            }
            k = k + 1;
        }//end    while(k<350)
        int bigSwaveLogic = (QRSPeaks[0] > mainWaveHeight / 3) && (QRSPeaks[1] > mainWaveHeight / 3) &&
                            (QRSPeaks[2] > mainWaveHeight / 3);
        if (((k > 281) && (QRSIso > localMaxPosition)) || bigSwaveLogic == 1 ||
            ((abs(QRSIso - localMaxPosition) < ISO_LENGTH_DIFF * 3) && (qrsIndex == 2) && (QRSPeaks[1] < 0))) {
            QRSIso = localMaxPosition - ISO_LENGTH_DIFF;
            // %  QRSIso=FIDMARK+floor((localMaxPosition-FIDMARK)/2);
        }
    }//end    if(QRSIso>280)//%%&&(pacedFlag>=BEAT_MS24)
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int QRSOffset = 0;
    if (((*SPeak != 0) && ((qrsIndex >= 2) || sWaveFlag == 1)) ||
        ((*SPeak > 198) && (*SPeak < QRSIso) && (*SPeak <= speakThrehold)))
        //%%  QRSOffset  = GetMinimmumAngle( beat,SPeak+BEAT_MS20,QRSIso );
    {
        if ((qrsIndex == 1) && (QRSPeaks[qrsIndex] < 0) && (QRSIso > BEAT_MS520)) {
            QRSIso = QRSIso - BEAT_MS20 * 2;
            QRSOffset = GetMinimmumAngleB(beat, *SPeak + BEAT_MS20 * 2, QRSIso);
        } else {
            // QRSOffset  = GetMinimmumAngle( beat,*SPeak+BEAT_MS20,QRSIso );
            if ((*SPeak > FIDMARK) && (qrsIndex > 1) && (fabsf(beat[(*SPeak)]) > 0.4 * fabsf(beat[FIDMARK]))) {
                QRSOffset = GetMinimmumAngleB(beat, *SPeak + BEAT_MS20, QRSIso);
            } else {
                QRSOffset = GetMinimmumAngleB(beat, *SPeak, QRSIso);
            }

        }
    } else {
        if (SSlope != 0)
            QRSOffset = GetMinimmumAngleB(beat, SSlope, QRSIso);
        //end
    }
// %%%%%%%%  big S Wave  %%%%%%%% rS morphology
    if ((*SPeak != 0) && (qrsIndex > 1) && (SSlope != 0) && (fabsf(beat[(int) (*SPeak)]) > 0.8 * fabsf(beat[FIDMARK])))
        QRSOffset = GetMinimmumAngleB(beat, SSlope, QRSIso);
    //end
    int QRSOffsetFirst = 0;
    if ((*SPeak > BEAT_MS480) && (fabsf(beat[(int) (*SPeak)]) > 2 * maxQRS)) {
        QRSOffsetFirst = GetMinimmumAngleB(beat, *SPeak + BEAT_MS20, QRSIso);
        if (QRSOffsetFirst > QRSOffset) {
            QRSOffset = QRSOffsetFirst;
        }
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// %%%%%%%%%%%%%%%%%%%%  withou S wave conditions  14:11 08/01 2019 ypz  %%%%%%%%%%%%%%%%%%%%
    if ((speakTrueFlag == 0) && ((SSlope != 0))) {
        if ((abs(QRSIsoSecond - localMaxPosition) < ISO_LENGTH_DIFF * 3) && (QRSPeaks[1] > 0) && (qrsIndex == 2) &&
            QRSPeaks[1] < 0)
            QRSIsoSecond = FIDMARK + (localMaxPosition - FIDMARK) / 2;
        // end
        if (sWaveFlag == 0)
            QRSOffset = GetMinimmumAngleB(beat, SSlope, QRSIsoSecond);
        // end
    }
//% %%%%%%%%%%%%%%%%%%%%%%%%%%% big and wide S wave  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  16:08 07/09 2017 ypz
    int index;
    int bigSwaveCount = 0;
    if ((*SPeak > FIDMARK) && (beat[(int) (*SPeak)] < 0) && (*SPeak < BEAT_MS500) &&
        (fabsf(beat[(int) (*SPeak)]) > mainWaveHeight / 2)) {
        index = *SPeak - BEAT_MS30;
        if (index < FIDMARK)
            index = *SPeak - BEAT_MS20;
        //end
        //    bigSwaveCount=0;
        while ((index > FIDMARK) && (index < (*SPeak + BEAT_MS20))) {
            if (fabsf(beat[index] - beat[(int) (*SPeak)]) < Crit)
                bigSwaveCount = bigSwaveCount + 1;
            //end
            index = index + 1;
        }
        if (bigSwaveCount >= BEAT_MS20)
            QRSOffset = GetMinimmumAngleB(beat, *SPeak + BEAT_MS20, QRSIso);
        //end
    }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ((qrsIndex == 1) && (QRSPeaks[0] > 0) && (beat[minPosition] < 0) && (QRSIso < minPosition) && (SSlope != 0) &&
        (minPosition > FIDMARK)) {
        int temPosition = (int) ((float) QRSIso + floorf((float) (minPosition - QRSIso) / 2));
        QRSOffset = GetMinimmumAngleB(beat, SSlope + BEAT_MS20 + 5, temPosition);
    }
    if ((*SPeak == 0) && (SpeakIndex < 280)) {
        *SPeak = SpeakIndex;
    }
    // end

    if (*SPeak != 0) {
        *rweight = 0;
        *speakheight = beat[(int) (*SPeak)];
        if (*SPeak != FIDMARK) {
            if (fabsf(beat[FIDMARK]) > fabsf(beat[(int) (*SPeak)])) {
                *rweight = beat[FIDMARK];
            } else {
                *rweight = beat[(int) (*SPeak)];
            }
        } else {
            *rweight = beat[(int) (*SPeak)];
        }
    }// end
    // %%%%%%%%%%%%% updated 18/03 2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  %  if(minSlopeRight> Crit/2)
    //  % %    QRSOffset  = GetMinimmumAngle( beat,SPeak+BEAT_MS20,QRSIso );
    //  %   end
    *onset = QRSstart;
    *offset = QRSOffset;
    if ((*offset == 0) || (*offset < FIDMARK))
        *offset = QRSfinish + 6;

    //end
    //%%  single QS morphology
    if ((fabsf(beat[FIDMARK]) >= mainWaveHeight * 0.9) && (beat[FIDMARK] < 0) && (QRSIso > BEAT_MS520))
        if (QRSOffset < (QRSfinish - 20))
            *offset = QRSfinish;
    //end
    //end
    // if (*offset < QRSPeaksPositions[qrsIndex - 1])
    //     *offset = QRSPeaksPositions[qrsIndex - 1];
    //end
    if ((*onset != 0) && (iosStart != 0)) {
        if (fabsf(beat[*onset]) < fabsf(beat[iosStart])) {
            *isoelectricLevel = pIsoLevelB(beat, iosStart);
        }
    }
    *rweight = mainWaveHeight;
    if ((*offset > BEAT_MS1000) && (QRSfinish < BEAT_MS1000)) {
        *offset = QRSfinish;
    }
// %%%%%%%%%%% Baseline Wander Correction  %%%%%%%%%%=
}


