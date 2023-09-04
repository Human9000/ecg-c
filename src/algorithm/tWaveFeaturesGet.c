#include "twavefeaturesget.h"
#include <math.h>

#define BEAT_MS20 10
#define BEAT_MS1000 500
#define BEAT_MS900 450
#define BEAT_MS960 480
#define BEAT_MS540 270
#define BEAT_MS400 200
#define BEAT_MS80 40
#define BEAT_MS10 5
#define limit 0.08
#define Fs BEAT_MS1000
#define FIDMARK 200
#define sampleRate BEAT_MS1000
#define dataLength sampleRate
#define M_PI 3.1415926535897

#pragma warning(disable:4244)



//TwaveFeaturesGet::TwaveFeaturesGet()
//{





//}

//TwaveFeaturesGet::~TwaveFeaturesGet()
//{




//}
//   //////////////////////////////////////////////////////////////////////
/// \briefvoid TwaveFeaturesGet::localMaxAlgorithm
/// \param float* array
/// \param int begin
/// \param int finish

void
localMaxAlgorithmT(float *array, int begin, int finish, float *tmax, int *maxPosition, float *tmin, int *minPosition) {
    // len=length(array);
    float D[dataLength] = {0};//zeros(1,len);
    float y[dataLength] = {0};//zeros(1,len);
    // float p[dataLength]={0};//zeros(1,len);
    *tmax = D[begin + 1];
    *tmin = D[begin + 1];
    *maxPosition = begin;
    *minPosition = begin;


    if (finish >= 500) {
        finish = 500;
    }

    if (begin >= 500) {
        begin = 500;
    }
    int ii = begin;
    while (ii < finish) {
        y[ii] = array[begin] + (float) (ii - begin) * (array[finish] - array[begin]) / (float) (finish - begin);
        D[ii] = fabsf(y[ii] - array[ii]);

        if (D[ii] > *tmax) {
            *tmax = D[ii];
            *maxPosition = ii;
        } else {
            if (D[ii] < *tmin) {
                *tmin = D[ii];
                *minPosition = ii;
            }
        }


        ii++;

    }


}
//   //////////////////////////////////////////////////////////////////////
/// \brief qrsFeatures::IsoCheckDiff
/// \param data
/// \param j
/// \param isoLength

int IsoCheckDiffT(float *data, int j, int isoLength, float tlimit) {

    // %   IsoCheck determines whether the amplitudes of a run
    // %   of data fall within a sufficiently small amplitude that
    // %	the run can be considered isoelectric.
    //   % BEAT_MS20=10;

    int i = j;
    int k = 0;
    // %%% search in the interval  from the biggest peak to 120ms
    while ((k < isoLength) && (fabsf(data[i + k] - data[i + k - 1]) < tlimit))
        k = k + 1;

    if ((k >= isoLength) && ((fabsf(data[i] - data[i + isoLength - 1])) < 3 * tlimit)) {
        return 0;
    } else {
        return 1;
    }
}

//   //////////////////////////////////////////////////////////////////////
// / \briefvoid  int   GetMinimmumAngle(float* beat,int leftmost,int rightmost )
// / \param float* beat
// / \param int leftmost
//   \int rightmost

int GetMinimmumAngleT(const float *beat, int leftmost, int rightmost) {

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
    int temPosition;
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
        if (leftX == 0) {
            // int  kk=0;
            // int  ll=kk;
            leftX = 1;
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
// / \briefvoid TwaveFeaturesGet::LowPassFilter
// / \param const float*X
// / \param float*Xpc
//
void LowPassFilterT(const float *X, float *Xpc) {
    // %  Linear integer lowpass filter
    unsigned int sample = 500;
    int Fpa = 50;
    int mpa = (int) (sample / Fpa);

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
        }
        data1[0] = X[ptr];

        ptr = ptr + 1;
    }

    int Tpa = (mpa - 1);
    int T = Tpa + 1;
    int j = 0;
    while (j < (dataLength - T)) {
        Xpc[j] = Xpc[j + T];
        j = j + 1;
    }

    j = dataLength - (T - 1);
    while (j < dataLength) {
        Xpc[j] = 0;
        j = j + 1;
    }

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int isoCheckT(const float *x, int istart, int len, float ISO_LIMIT) {
    float res = 0;
    // %tmp = x(istart:istart+len-1);
    float temp[10] = {0};
    // temp1= ecgCopy(x,istart,istart+len-1);
    float tempMin = x[istart];
    float tempMax = x[istart];
    int ii = istart;
    while (ii < len) {
        temp[ii] = x[ii];
        if (x[ii] > tempMax) {
            tempMax = x[ii];
        } else if ((x[ii] < tempMin)) {
            tempMin = x[ii];
        }

        ii++;
    }
    //%t1=max(tmp);
    //    tt1= ECGMax(temp,len);
    // %t2=min(tmp);
    //     tt2= ECGMin(temp,len );
    // %if max(tmp) - min(tmp) > ISO_LIMIT
    if ((tempMax - tempMin) > ISO_LIMIT) {
        res = 0;
    } else if ((fabsf(tempMax) < ISO_LIMIT)) {
        res = 1;
    }
    return (int) res;
}

//   //////////////////////////////////////////////////////////////////////
// / \brief void TwaveFeaturesGet::tWavesFeaturesGet
// / \param int onset
// / \param int rr
//
// / \param int offset
// / \param int QPeak
// / \inint SPeak

void tWavesFeaturesGet(float *beat, int onset, int rr, int offset, int SPeak, float isoelectricLevel, float tLimit,
                       float Crit,
                       twav_data_t *tPoint)//int*tFirstPeak,int*tSecondPeak,int*type,float*tHeight,int*Tonset,int*tOffset,float*tSum)
{

    tPoint->tFirstPeak = 0;
    tPoint->tSecondPeak = 0;
    tPoint->type = 0;
    tPoint->tHeight = 0;
    tPoint->Tonset = 0;
    tPoint->tOffset = 0;
    tPoint->tAera = 0;
    tPoint->qtInterval = 0;


    int ii = 0;

    float Xpc[sampleRate] = {0};
    LowPassFilterT(beat, Xpc);

    // %%%%%%%%%%%%%%%%%%%%%%% Wings Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    double RRav = (double) rr / BEAT_MS1000;
    float HR = (float) (60 / RRav);
    int SToffset;
    if (HR >= 140) {
        SToffset = (int) (offset + floor(0.04 * Fs));  //%% 40ms
    } else if (HR >= 120) {
        SToffset = (int) (offset + round(0.06 * Fs));
    } else if (HR >= 110) {
        SToffset = (int) (offset + round(0.064 * Fs));
    } else if (HR >= 100) {
        SToffset = (int) (offset + round(0.072 * Fs));
    } else if (HR < 100) {
        SToffset = (int) (offset + round(0.08 * Fs));
    }// end

    int tInitial = 20;
    int tBoundary = 500;
    int tWindowEnd = (int) round(0.4 * sqrt(RRav) * Fs * 8 / 7 + onset) + 2 * BEAT_MS20;
    //   int Tend=onset+round(0.55*(sqrt(RRav)*Fs)); // %% National University of defence and technology ......
    if (tLimit > BEAT_MS900)
        tLimit = BEAT_MS900;
    if ((tLimit != 0) && (rr < BEAT_MS540) && (tLimit < BEAT_MS900) && (SPeak != 0))
        tWindowEnd = (int) (tLimit - BEAT_MS20);

    if (tWindowEnd >= tBoundary)
        tWindowEnd = tBoundary - tInitial * 2;


    ii = offset + tInitial / 4;
    if (ii > tWindowEnd)
        ii = (int) floorf(FIDMARK + (float) (tWindowEnd - FIDMARK) / 2);

    float tWing[sampleRate] = {0};


    //  $$$$$$$$$$$$$$$$$$$$$$$$$
    if ((tWindowEnd < BEAT_MS1000) && ((onset + rr) <= BEAT_MS1000))

        //if(Ponset>0){
        //  tWindowEnd=BEAT_MS400+rr;
        // }
        // else{
        tWindowEnd = onset + rr;
        // }
        ////end
    else {
        ii = FIDMARK + BEAT_MS80;
        float maxQRS = (Xpc[ii]);
        float minQRS = (Xpc[ii]);
        int maxPosition = ii;
        int minqrsPosition = ii;
        // looking for local max peak
        while (ii < FIDMARK + BEAT_MS400) {
            if (((Xpc[ii]) > maxQRS)) {
                maxQRS = fabsf(Xpc[ii]);
                maxPosition = ii;
            } else {
                if (((Xpc[ii]) < minQRS))//&&((fabs(beat[ii])>beat[ii+BEAT_MS10])&&(fabs(beat[ii])>beat[ii-BEAT_MS10])))
                {
                    minQRS = fabsf(Xpc[ii]);
                    minqrsPosition = ii;
                }
            }//end
            ii = ii + 1;
        }//end

        if (minqrsPosition > maxPosition) {
            ii = minqrsPosition;
        } else {
            ii = maxPosition;
        }

        // looking for ISOelectric
        int miniPosition = 0;
        while ((ii < BEAT_MS960) && (isoCheckT(Xpc, ii, (int) (0.02 * Fs), Crit / 2) == 0)) {
            ii = ii + 1;
        }
        if (ii < BEAT_MS960)
            miniPosition = (int) (ii + 0.02 * Fs);
        else {
            ii = maxPosition;//% FIDMARK+0.16*fs;
            while ((ii < BEAT_MS960) && (isoCheckT(Xpc, ii, (int) (0.02 * Fs), Crit) == 0)) {
                ii = ii + 1;
            }
            if (ii < BEAT_MS960)
                miniPosition = (int) (ii + 0.02 * Fs);
            else {
                ii = maxPosition;//% FIDMARK+0.16*fs;
                while ((ii < BEAT_MS960) && (isoCheckT(Xpc, ii, (int) (0.02 * Fs), Crit * 2) == 0)) {
                    ii = ii + 1;
                }//end
                if (ii < BEAT_MS960)
                    miniPosition = (int) (ii + 0.02 * Fs);
                else
                    miniPosition = maxPosition + BEAT_MS80;
            }//end
        }//end

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tWindowEnd = miniPosition;
    }// end
    if (tWindowEnd >= 460)
        tWindowEnd = 460;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float tWinMax = 0;
    int tWinMaxPosi = 0;
    float tWinMin = 0;
    int tWinMinPosi = 0;

    ii = offset + tInitial / 4;
    float wingLeft;
    float wingRight;
    while ((ii < tWindowEnd) && (ii < (tBoundary - tInitial))) {
        if (Xpc[ii] > tWinMax) {
            tWinMax = Xpc[ii];
            tWinMaxPosi = ii;
        } else if (Xpc[ii] < tWinMin) {
            tWinMin = Xpc[ii];
            tWinMinPosi = ii;
        }

        int tempii = ii - tInitial;
        if (tempii < 0)
            tempii = 0;
        wingLeft = Xpc[tempii] - Xpc[ii];
        wingRight = Xpc[ii] - Xpc[ii + tInitial];
        tWing[ii] = wingLeft * wingRight;
        ii++;
    }

    int tInterval = 10;
    float biphasic = 0.8f;
    if (offset < 0)
        offset = 0;

    int tWingStart;
    ii = offset + tInitial * 2;
    if (Xpc[offset] < 0)     //%% S wave offset should not be high over zero line  15:16 11/02 2018
        tWingStart = offset + tInitial * 2;
    else {
        tWingStart = offset;
        ii = offset;
    }
    float tsmallest = tWing[ii];
    int smallest = ii;
    float tsmaller = tWing[ii];
    int smaller = ii;
    //  int tStart=ii;

    int tlag = 5;
    ii = ii + 2;
    while (ii < tWindowEnd) // 16:01 11/02 2018 ypz
    {
        if ((tWing[ii] < tsmallest) && ((tWing[ii] - tWing[ii - 2]) * (tWing[ii + 2] - tWing[ii]) < 0) &&
            (tWing[ii] < 0)) {
            tsmaller = tsmallest;
            smaller = smallest;
            tsmallest = tWing[ii];
            smallest = ii;
            ii = ii + tlag * 2;
        } else {
            if ((tWing[ii] < tsmaller) && ((tWing[ii] - tWing[ii - 2]) * (tWing[ii + 2] - tWing[ii]) < 0) &&
                (tWing[ii] < 0)) {
                tsmaller = tWing[ii];
                smaller = ii;
                ii = ii + tlag;
            }
        }
        ii = ii + 1;
    }


    tPoint->Tonset = 0;
    if (smaller > tWingStart + tlag * 2)//%% greater than initial value
    {
        if ((smaller > smallest) && (fabsf(Xpc[smaller]) > fabsf(Xpc[smallest]))) {
            float tempqq = fabsf(Xpc[smallest] / Xpc[smaller]);
            if ((fabsf(Xpc[smallest] / Xpc[smaller]) > biphasic)) {
                tPoint->tFirstPeak = smaller;
                tPoint->tSecondPeak = smallest;
            } else {
                tPoint->Tonset = smallest;
                tPoint->tFirstPeak = smaller;
            }
        } else {
            if ((smaller > smallest) && (fabsf(Xpc[smaller]) < fabsf(Xpc[smallest]))) {
                if ((fabsf(Xpc[smaller] / Xpc[smallest]) > biphasic)) {
                    tPoint->tFirstPeak = smallest;
                    tPoint->tSecondPeak = smaller;
                } else {
                    tPoint->tOffset = smaller;
                    tPoint->tFirstPeak = smallest;
                }
            } else {
                if ((smaller < smallest) && (fabsf(Xpc[smallest]) > fabsf(Xpc[smaller]))) {
                    if ((fabsf(Xpc[smaller] / Xpc[smallest]) > biphasic)) {
                        tPoint->tFirstPeak = smallest;
                        tPoint->tSecondPeak = smaller;
                    } else {
                        tPoint->Tonset = smaller;
                        tPoint->tFirstPeak = smallest;
                    }
                } else {

                    if ((smaller < smallest) && (fabsf(Xpc[smallest]) < fabsf(Xpc[smaller]))) {
                        if ((fabsf(Xpc[smallest] / Xpc[smaller]) > biphasic)) {
                            tPoint->tFirstPeak = smaller;
                            tPoint->tSecondPeak = smallest;
                        } else {
                            tPoint->tFirstPeak = smaller;
                            tPoint->tSecondPeak = 0;
                        }
                    }

                }

            }

        }

    } else {
        smaller = 0;
        if (smallest > tWingStart)
            tPoint->tFirstPeak = smallest;
    }

    if ((tPoint->tFirstPeak == 0) && (tPoint->tSecondPeak == 0))
        tPoint->tFirstPeak = tWindowEnd;
    //end

    int Toffset = 0;
    int QRSOffset;
    if ((tPoint->tFirstPeak != 0)) {   //for level
        // ii=tPoint->tFirstPeak;

        if (smaller > smallest)
            ii = smaller;
        else
            ii = smallest;
        int centerPoint = ii;

        while ((ii < tWindowEnd) && (fabsf(Xpc[centerPoint] - Xpc[ii + 1]) < Crit * 1.5))
            ii++;
        if (ii < tWindowEnd) {
            QRSOffset = GetMinimmumAngleT(Xpc, ii, tWindowEnd - tInterval / 2);
        } else {
            QRSOffset = GetMinimmumAngleT(Xpc, (int) round(tPoint->tFirstPeak + tInterval), tWindowEnd - tInterval);
        }


        float tmax = 0;
        int tmaxPosition = 0;
        float tmin = 0;
        int tminPosition = 0;
        if (QRSOffset != 0) {
            Toffset = QRSOffset;
        } else {
            tBoundary = (int) (tWindowEnd + floor((BEAT_MS1000 - tWindowEnd) * 0.8));
            if (((float) tBoundary < tLimit) && tBoundary < BEAT_MS1000)//(tLimit<BEAT_MS1000)
                tWindowEnd = tBoundary;
            //end
            localMaxAlgorithmT(Xpc, (int) round(tPoint->tFirstPeak + BEAT_MS10), (int) round(tWindowEnd), &tmax,
                               &tmaxPosition,
                               &tmin, &tminPosition);
            Toffset = tmaxPosition;
        }

        if (tPoint->Tonset == 0) {
            if ((tPoint->tSecondPeak != 0) && (tPoint->tSecondPeak < tPoint->tFirstPeak) &&
                (tPoint->tSecondPeak > offset)) {
                QRSOffset = GetMinimmumAngleT(Xpc, offset, tPoint->tSecondPeak);
                tPoint->Tonset = QRSOffset;
            } else if ((tPoint->tSecondPeak == 0) && (smaller > 0) && (smaller < tPoint->tFirstPeak))
                tPoint->Tonset = smaller;
            //end
            //end

            if ((tPoint->tFirstPeak > 0) && (offset > 0) && (tPoint->tFirstPeak > offset) && (tPoint->Tonset == 0))
                // [ tmax,tmaxPosition,tmin,tminPosition ] = localMaxAlgorithm( ecgdata,round(offset+tlag*4),round(tFirstPeak) );%% 40ms after QRS j point
                localMaxAlgorithmT(Xpc, (int) round(offset + tlag * 4), (int) round(tPoint->tFirstPeak), &tmax,
                                   &tmaxPosition,
                                   &tmin, &tminPosition);
            tPoint->Tonset = tmaxPosition;
            //end

            if (SToffset > tPoint->Tonset)
                tPoint->Tonset = SToffset;    //    15:36 30/01 2018 ypz Mao Ling ,National University of Defence Technology
            //end
        }

    }// end


    // tSum=0;
    if (tPoint->tFirstPeak > 15) {
        int k = tPoint->tFirstPeak - 15;
        float sum = 0;
        while (k <= (tPoint->tFirstPeak + 20)) {
            if (beat[k] > 0) {
                sum = sum + beat[k];
            } else {
                sum = sum + fabsf(beat[k]);
            }
            k++;
        }
        tPoint->tAera = sum;
    }


    tPoint->tOffset = Toffset;
    //  %  type=result.type;
    if (tPoint->tFirstPeak < (tWindowEnd)) {
        if ((tPoint->tFirstPeak != 0) && (Toffset != 0)) {
            tPoint->tHeight = beat[tPoint->tFirstPeak];
            tPoint->tOffset = Toffset;
        }
    } else if (tWinMinPosi > 0) {
        tPoint->tHeight = beat[tWinMinPosi];
        tPoint->tOffset = Toffset;
        tPoint->tFirstPeak = tWinMinPosi;
    }


    if (tPoint->tHeight > 0) {
        tPoint->type = 1;
    } else {
        tPoint->type = 0;
    }

    if ((tPoint->Tonset != 0) && (tPoint->tOffset != 0)) {
        tPoint->qtInterval = tPoint->tOffset - onset;
    }


}
