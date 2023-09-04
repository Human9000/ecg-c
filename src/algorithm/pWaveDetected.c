#include "pWaveDetected.h"
#include <math.h>
#include "utility.h"

#define BEAT_MS16 8
#define BEAT_MS20 10
#define BEAT_MS40 20
#define BEAT_MS110 55
#define BEAT_MS120 60
#define BEAT_MS130 65
#define BEAT_MS150 75
#define BEAT_MS270 135
#define BEAT_MS220 110
#define BEAT_MS200 100
#define BEAT_MS210 105
#define BEAT_MS400 200
#define BEAT_MS600 300
#define maxInteger 32767
#define minInteger (-32768)

#define EPSILON 0.000001
#define negaThreshold (1/15)
#define pLimit 200
#define pError 10
#define limit 0.085
#define segmentLen 50

#define pdataLength 260
#define level 5
#define Fs 500
#define FIDMARK 200
#define sampleRate 500
#define modulusDist 75
#define pLength pLimit
#define BEAT_MS140 70

#define sampleRate 500
#define dataLength sampleRate
#define len 200



//pWaveDetected::pWaveDetected()
//{
//  int jj=0;
//  int ii=0;

//   result.peak1=0;
//   result.type=type_0;
//   result.peak2=0;
//   result.onset=0;
//   result.offset=0;
//   result.weight=0;

//  ii=0;
//  while(ii<6)
//  {
//    Presult[ii].onset=0;
//    Presult[ii].offset=0;
//    Presult[ii].type=type_0;
//    Presult[ii].coeff=0;
//    Presult[ii].Peak=0;
//    Presult[ii].weight=0;
//     ii++;
//  }


//  ii=0;

//  while(ii<5)
//  {
//     jj=0;
//     while(jj<500)
//      {
//       swa[ii][jj]=0;
//       swd[ii][jj]=0;
//       jj++;

//      }
//      ii++;
//    }

//}

//pWaveDetected::~pWaveDetected()
//{



//}

//   //////////////////////////////////////////////////////////////////////
// / \brief void pWaveDetected:: pWaveRoutine
// / \param int FIDMARK
// / \param int SF
// / \param int*pPonset
// / \param int*pPoffsetn
// / \int*Ppos

int pWaveIsoCheck(const float *x, int istart, int pplen, float ISO_LIMIT) {
    // utility pTemp;
    float temp[25] = {0};
    int ii = 0;
    int jj = istart;
    while (jj < (istart + pplen - 1)) {
        temp[ii] = x[istart + jj];
        ii++;
        jj++;
    }
    if ((ECGMax(temp, pplen) - ECGMin(temp, pplen)) > ISO_LIMIT)
        return 0;
    else
        return 1;

}


//   //////////////////////////////////////////////////////////////////////
// / \brief void pWaveDetected:: pWaveRoutine
// / \param int FIDMARK
// / \param int SF
// / \param int*pPonset
// / \param int*pPoffsetn
// / \int*Ppos

void pWaveRoutine(float *Xpc, int *Ponset, int *Poffset, int *Ppos) {
    int findmark = FIDMARK;
    int fs = sampleRate;
    float dx[FIDMARK] = {0};
    // utility pTemp;
    ecgDiff(Xpc, dx, BEAT_MS400);
    int ii = FIDMARK - 0.25 * fs;
    if (ii < 0) {
        ii = 0;
    }
    float maxflue = 0;
    int maxflueI = ii;
    float temp[25] = {0};
    while (ii < findmark - 0.1 * fs) {
        ecgCopy(dx, temp, ii, (int) (ii + 0.05 * fs));
        if (maxflue < (ECGMax(temp, (int) (0.05 * fs)) - ECGMin(temp, (int) (0.05 * fs)))) {
            maxflue = ECGMax(temp, (int) (0.05 * fs)) - ECGMin(temp, (int) (0.05 * fs));
            maxflueI = ii;
        }
        ii = ii + 1;
    }


    // % ??????,?????????   pWaveIsoCheck(float*x,int istart,int pplen,float ISO_LIMIT)
    ii = maxflueI;
    float x[25] = {0};
    ecgCopy(Xpc, x, ii, (int) (ii + 0.05 * fs));
    maxflue = ECGMax(x, (int) (0.05 * fs)) - ECGMin(x, (int) (0.05 * fs));
    while ((ii < (maxflueI + 0.05 * fs)) && (pWaveIsoCheck(Xpc, ii, (int) (0.05 * fs), maxflue / 2) == 0)) {
        ii = ii + 1;
    }
    *Poffset = (int) (ii + 0.05 * fs);


    // % ???????????,???????  pWaveIsoCheck(Xpc,ii,0.05*fs,maxflue/2)
    if (ii == maxflueI + 0.05 * fs) {
        ii = maxflueI;
        while (ii < maxflueI + 0.05 * fs && pWaveIsoCheck(Xpc, ii, (int) (0.03 * fs), maxflue / 2) == 0) {
            ii = ii + 1;
        }
        *Poffset = (int) (ii + 0.03 * fs);
    }

    //% ????????,
    ii = maxflueI;
    while (ii > (maxflueI - 0.1 * fs) && (ii > 0) && (pWaveIsoCheck(Xpc, ii, (int) (0.05 * fs), maxflue / 2) == 0)) {
        ii = ii - 1;
    }


    *Ponset = ii;
    //%???????
    if (ii == (maxflueI - 0.1 * fs)) {
        ii = maxflueI;
        while ((ii > maxflueI - 0.1 * fs) && ii > 0 && (pWaveIsoCheck(Xpc, ii, (int) (0.03 * fs), maxflue / 2) == 0)) {
            ii = ii - 1;
        }
        *Ponset = ii;
    }

    if (*Ponset < 0)
        *Ponset = 0;

    float z = (Xpc[*Ponset] + Xpc[*Poffset]) / 2;
    ii = *Ponset;
    int bb = 0;
    float ecgTempMax = fabsf(Xpc[ii] - z);
    while (ii < *Poffset) {
        if (fabsf(Xpc[ii] - z) > ecgTempMax) {
            ecgTempMax = fabsf(Xpc[ii] - z);
            bb = ii;
        }
        ii++;
    }

    *Ppos = bb - 1;

}


//   //////////////////////////////////////////////////////////////////////
// / \brief void pWaveDetected::ECGMaxPosition
// / \param const float *sourceBeat
// / \param int len
// / \param float *returnMax
// / \param int *position

int threeDataMax(float pWeight1, float pWeight2, float pWeight3) {
    float pMax;
    if ((pWeight1 > pWeight2 && pWeight1 > pWeight3)) {
        pMax = 2;
    } else {
        if (pWeight3 > pWeight1 && pWeight3 > pWeight2) {
            pMax = 4;
        } else {
            pMax = 3;
        }
    }

    return (int) pMax;
}

//   //////////////////////////////////////////////////////////////////////
// / \brief void pWaveDetected::ECGMaxPosition
// / \param const float *sourceBeat
// / \param int len
// / \param float *returnMax
// / \param int *position

void ECGMaxPosition(const float *sourceBeat, int length, float *returnMax, int *position) {
    int tempPosition;
    if (length >= 1) {
        *returnMax = sourceBeat[0];
        int kk = 1;
        tempPosition = kk;
        while (kk < length) {
            if (sourceBeat[kk] > *returnMax) {
                *returnMax = sourceBeat[kk];
                tempPosition = kk;
            }
            kk = kk + 1;
        }
    } else {
        *returnMax = 0;
        *position = 0;
        return;
    }
    *position = tempPosition;

}


//   //////////////////////////////////////////////////////////////////////
/// \brief  void pWaveDetected::bufferCopy
/// \param int start
/// \param int finish
/// \param float*returnBuffer
/// \param

void bufferCopy(const float *sourceBeat, int start, int finish, float *returnBuffer) {

    int kk = start;
    int mm = 0;
    while (kk < finish) {
        returnBuffer[mm] = sourceBeat[kk];
        mm = mm + 1;
        kk = kk + 1;

    }

}

//   //////////////////////////////////////////////////////////////////////
/// \brief  int pWaveDetected::localMaxAlgorithm
/// \param  begin
/// \param  finish
/// \param

int localMaxAlgorithm(float *array, int begin, int finish) {

    float D[200] = {0};
    float y[200] = {0};
    float max = -32768;
    float min = 32767;
    int maxPosition = begin;
    int minPosition = begin;

    int i = begin;
    while (i < finish) {
        y[i] = array[begin] + (float) (i - begin) * (array[finish] - array[begin]) / (float) (finish - begin);
        D[i] = fabsf(y[i] - array[i]);

        if (D[i] > max) {
            max = D[i];
            maxPosition = i;
        } else {
            if (D[i] < min) {
                min = D[i];
                minPosition = i;
            }
        }

        i = i + 1;
    }//end   while(i<finish

    return maxPosition;
}

//   //////////////////////////////////////////////////////////////////////
/// \brief  pWaveDetected::pSumCalculate
/// \param  ecgdata
/// \param  Ponset
/// \param  Poffset
/// \param  Psum

void pSumCalculate(float *ecgdata, int Ponset, int Poffset, float *Psum) {
    *Psum = 0;
    int i = 0;
    if (Ponset != 0) {
        i = Ponset;
        while (i < Poffset) {
            if (ecgdata[i] > 0)
                *Psum += ecgdata[i];
            else
                *Psum += fabsf(ecgdata[i]);
            i = i + 1;
        }
    }
}

//   //////////////////////////////////////////////////////////////////////
/// \brief  pWaveDetected::localMaxAlgorithm
/// \param array
/// \param begin
/// \param finish
/// \param max
/// \param  maxPosition
/// \param  min
/// \param minPosition


void localMaxAlgorithm1(float *array, int begin, int finish, float max, int maxPosition, float min, int minPosition) {

    float D[len];
    float y[len];
    float p[len];

    int ii = 0;
    while (ii < len) {
        D[ii] = 0;
        y[ii] = 0;
        p[ii] = 0;
        ii++;
    }

    max = D[begin + 1];
    min = D[begin + 1];
    maxPosition = begin;
    minPosition = begin;

    ii = begin;
    while (ii <= finish) {
        y[ii] = array[begin] + (float) (ii - begin) * (array[finish] - array[begin]) / (float) (finish - begin);
        D[ii] = fabsf(y[ii] - array[ii]);

        if (D[ii] > max) {
            max = D[ii];
            maxPosition = ii;
        } else if (D[ii] < min) {
            min = D[ii];
            minPosition = ii;
        }

        ii = ii + 1;
    }

}

//   //////////////////////////////////////////////////////////////////////
/// \briefpWaveDetected::LowPassFilter
/// \param  X
/// \param  sample
/// \param  Xpc


void LowPassFilterP(const float *X, float *Xpc) {
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

//   //////////////////////////////////////////////////////////////////////
/// \brief pWaveDetected:: PWaveletDetect
/// \param beat
/// \param onset
/// \param rr
/// \param isoelectricLevel
/// \param  preToffset
/// \param  Ponset
/// \param pPeak
/// \param Poffset
///

void PWaveletDetect(const float *beat, unsigned int onset, int rr, float isoelectricLevel, int *Ponset, float *pPeak,
                    int *Poffset, int *pPeakPositionFirst, int *pPeakPositionSecond, int *type, float *Psum, int init) {

    pwav_data_t2 result;
    pwav_data_t Presult[6];
    float swa[5][(pdataLength + 5)];
    float swd[5][(pdataLength + 5)];

    // if(init)
    //  {
    //  int jj=0;
    unsigned ii = 0;

    result.peak1 = 0;
    result.type = 0;//type_0;
    result.peak2 = 0;
    result.onset = 0;
    result.offset = 0;
    result.weight = 0;

    ii = 0;
    while (ii < 6) {
        Presult[ii].onset = 0;
        Presult[ii].offset = 0;
        Presult[ii].type = 0;//type_0;
        Presult[ii].coeff = 0;
        Presult[ii].Peak = 0;
        Presult[ii].weight = 0;
        ii++;
    }


    ii = 0;

    while (ii < 5) {
        int jj = 0;
        while (jj < (pdataLength + 5)) {
            swa[ii][jj] = 0;
            swd[ii][jj] = 0;
            jj++;

        }
        ii++;
    }

    // }



    *pPeak = 0;
    *pPeakPositionFirst = 0;
    *pPeakPositionSecond = 0;
    *type = 0;
    *Psum = 0;

    float Xpc[sampleRate];
    // utility pTemp;
    LowPassFilter(beat, Xpc);

    // % wavelet cofficients
    // % low pass 1/4 3/4 3/4 1/4
    // % high pass -1/4 -3/4 3/4 1/4
    // %
    // % quadratic B spline wavelet
    // % Pairs of Maximum moduli

    unsigned int j, jj, jj0, jj1, jj2, jj3;
    unsigned int i;


    for (ii = 0; ii < (pdataLength - 2); ii++) {
        jj0 = ii + 3 - 0;//-qPow(2,0)*0;
        jj1 = ii + 3 - 1;//-qPow(2,0)*1;
        jj2 = ii + 3 - 2;//-qPow(2,0)*2;
        jj3 = ii + 3 - 3;//-qPow(2,0)*3;

        swa[0][ii + 3] = 0.25f * Xpc[jj0] + 0.75f * Xpc[jj1] + 0.75f * Xpc[jj2] + 0.25f * Xpc[jj3];
        swd[0][ii + 3] = -0.25f * Xpc[jj0] - 0.75f * Xpc[jj1] + 0.75f * Xpc[jj2] + 0.25f * Xpc[jj3];
    }

    j = 1;
    while (j < level) {
        int pPow = 1;
        unsigned int k = 0;
        while (k < j) {
            pPow = pPow * 2;
            k++;
        }

        for (i = 1; i < pdataLength - 47; i++) {

            jj0 = i + 47 - pPow * 0;//qPow(2,j)*0;
            jj1 = i + 47 - pPow * 1;//qPow(2,j)*1;
            jj2 = i + 47 - pPow * 2;//qPow(2,j)*2;
            jj3 = i + 47 - pPow * 3;//qPow(2,j)*3;

            swa[j][i + 47] = 0.25f * swa[j - 1][jj0] + 0.75f * swa[j - 1][jj1] + 0.75f * swa[j - 1][jj2] +
                             0.25f * swa[j - 1][jj3];
            swd[j][i + 47] = -0.25f * swa[j - 1][jj0] - 0.75f * swa[j - 1][jj1] + 0.75f * swa[j - 1][jj2] +
                             0.25f * swa[j - 1][jj3];
        }

        j = j + 1;

    }

    //negative minimum-positive maximum pair(modulus maximum pair (MMP))

    float pPositive[level][pdataLength] = {{0}};

    j = 0;
    while (j < level) {
        ii = 1;
        while (ii < pdataLength) {
            if (swd[j][ii] > 0) {
                pPositive[j][ii] = swd[j][ii];
            } else {
                pPositive[j][ii] = 0;
            }
            ii = ii + 1;
        }
        j = j + 1;
    }

//%%%%%%%%%%%%%%%% first derivate  %%%%%%%%%%%%%%%%%%%%%%%%

    float pPositiveGreat[level][pdataLength - 1] = {{0}};
    ii = 0;
    while (ii < level) {
        j = 0;
        while (j < (pdataLength - 1)) {
            if ((pPositive[ii][j] - pPositive[ii][j + 1]) < 0) {
                pPositiveGreat[ii][j] = 1;
            } else {
                pPositiveGreat[ii][j] = 0;
            }
            j++;
        }
        ii++;
    }

//pddw(:,2:points-1)=((pdw(:,1:dataLength-2)-pdw(:,2:dataLength-1))>0);

    float pPosiGreatest[level][pdataLength] = {{0}};
    ii = 0;
    while (ii < level) {
        j = 0;
        while (j < (pdataLength - 2)) {
            if ((pPositiveGreat[ii][j] - pPositiveGreat[ii][j + 1]) > 0) {
                pPosiGreatest[ii][j + 1] = 1;
            } else {
                pPosiGreatest[ii][j + 1] = 0;
            }
            j++;
        }
        ii++;
    }

// negw=swd.*(swd<0);

    float pNegative[level][pdataLength] = {{0}};
    ii = 0;
    while (ii < level) {
        j = 0;
        while (j < pdataLength) {
            if (swd[ii][j] < 0) {
                pNegative[ii][j] = swd[ii][j];
            } else {
                pNegative[ii][j] = 0;
            }
            j++;
        }
        ii++;
    }

//ndw=((negw(:,1:points-1)-negw(:,2:points))>0)

    float pNegativeGreat[level][pdataLength - 1] = {{0}};

    ii = 0;
    while (ii < level) {
        jj = 0;
        while (jj < (pdataLength - 1)) {
            if ((pNegative[ii][jj] - pNegative[ii][jj + 1]) > 0) {
                pNegativeGreat[ii][jj] = 1;
            } else {
                pNegativeGreat[ii][jj] = 0;
            }
            jj++;
        }
        ii++;
    }

    //   %?????
    //   % nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);

    float pNegaGreatest[level][pdataLength - 1] = {{0}};
    ii = 0;
    while (ii < level) {
        jj = 0;
        while (jj < (pdataLength - 2)) {
            if ((pNegativeGreat[ii][jj] - pNegativeGreat[ii][jj + 1]) > 0) {
                pNegaGreatest[ii][jj + 1] = 1;
            } else {
                pNegaGreatest[ii][jj + 1] = 0;
            }
            jj++;
        }
        ii++;
    }
//      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//      %???
//      %ddw=pddw|nddw;
//      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float pPosiNegative[level][pdataLength] = {{0}};
    ii = 0;
    while (ii < level) {
        jj = 0;
        while (jj < pdataLength) {
            if (pNegaGreatest[ii][jj] > EPSILON || pPosiGreatest[ii][jj] > EPSILON) {
                pPosiNegative[ii][jj] = 1;
            } else {
                pPosiNegative[ii][jj] = 0;
            }
            jj++;
        }
        ii++;
    }

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        % ddw(:,1)=1;
//        % ddw(:,points)=1;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 0;
    while (ii < level) {
        pPosiNegative[ii][0] = 1;
        pPosiNegative[ii][pdataLength - 1] = 1;
        ii++;
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %???????,????0
// %  wpeak=ddw.*swd;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float pWaveletPeak[level][pdataLength] = {{0}};
    ii = 0;
    while (ii < level) {
        jj = 0;
        while (jj < pdataLength) {
            pWaveletPeak[ii][jj] = pPosiNegative[ii][jj] * swd[ii][jj];
            jj++;
        }
        ii++;
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %  wpeak(:,1)=wpeak(:,1)+1e-10;
// %  wpeak(:,points)=wpeak(:,points)+1e-10;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ii = 0;
    while (ii < level) {
        pWaveletPeak[ii][0] = pWaveletPeak[ii][0] + 1e-10f;
        pWaveletPeak[ii][pdataLength - 1] = pWaveletPeak[ii][pdataLength - 1] + 1e-10f;
        ii++;//kk=kk+1;
    }//end

// % Mj1=wpeak(1,:);
// % Mj2=wpeak(2,:);
// % Mj3=wpeak(3,:);
//% %Mj4=wpeak(4,:);
// % Mj4=Mj3;


    float pMj4[pdataLength] = {0};//zeros(1,dataLength);
    ii = 0;
    while (ii < pdataLength) {
        //% pMj4(kk)=pWaveletPeak(4,kk);
        pMj4[ii] = pWaveletPeak[3][ii];//%%%%12:46 13/09 2018
        ii++;
    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%5
// %Mj5=wpeak(5,:);
// %scale 3 minimun-maximum modulus pairs (MMP)
// % Mj4(pLength)=0;
// % Mj4(pLength-1)=0;

    pMj4[pLength - 1] = 0;
    pMj4[pLength - 2] = 0;


    // % posi=Mj4.*(Mj4>0);
    // % posi=posi(1:pwaveEnd);
    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float pWaveletPosi[pLength] = {0};//zeros(1,pLength);
    ii = 0;
    while (ii < pLength) {
        if (pMj4[ii] > 0) {
            pWaveletPosi[ii] = pMj4[ii];
        }
        ii++;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %positive maximum average
    //  points=pLength;%%
    // % thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))
    // +max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float array[50] = {0};//zeros(1,round(points/4));
    ii = 0;
    jj = 0;
    while (ii < 50) {
        array[jj] = pWaveletPosi[ii];
        ii++;
        jj++;
    }

    float temp1 = ECGMax(array, 50);

    ii = 0;
    jj = 50;
    while (jj < 100) {
        array[ii] = pWaveletPosi[jj];
        ii++;
        jj++;
    }
    float temp2 = ECGMax(array, 50);

    ii = 0;
    jj = 100;
    while (jj < 150) {
        array[ii] = pWaveletPosi[jj];
        ii++;
        jj++;
    }
    float temp3 = ECGMax(array, 50);

    ii = 0;
    jj = 150;
    while (jj < 200) {
        array[ii] = pWaveletPosi[jj];
        jj++;
        ii++;
    }
    float temp4 = ECGMax(array, 50);
    float temp = (temp1 + temp2 + temp3 + temp4) / 4;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// % posi=(posi>thposi*1/10);%%  1/3
    float waveletPosi[pLength] = {0};
    ii = 0;
    while (ii < pLimit) {
        if (pWaveletPosi[ii] > (temp * 0.1)) {
            waveletPosi[ii] = 1;
        } else {
            waveletPosi[ii] = 0;
        }
        ii++;
    }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//      % nega=Mj4.*(Mj4<0);
//      % nega=nega(1:pwaveEnd);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float pWaveletNega[pLength] = {0};//zeros(1,pLength);
    ii = 0;
    while (ii < pLength) {
        if (pMj4[ii] < 0) {
            pWaveletNega[ii] = pMj4[ii];
        }
        ii++;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %negative maximum average
    // %thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))
    // +min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round((points)/4):4*round(points/4))))/4;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ii = 0;
    jj = 0;
    while (jj < segmentLen) {
        array[ii] = pWaveletNega[jj];
        ii++;
        jj++;
    }
    temp1 = ECGMin(array, segmentLen);

    ii = 0;
    jj = segmentLen;
    while (jj < (2 * segmentLen)) {
        array[ii] = pWaveletNega[jj];
        ii++;
        jj++;
    }
    temp2 = ECGMin(array, segmentLen);

    ii = 0;
    jj = 2 * segmentLen;//2*round(points/4);
    while (jj < 3 * segmentLen) {
        array[ii] = pWaveletNega[jj];
        ii++;
        jj++;
    }
    temp3 = ECGMin(array, 3 * segmentLen);

    ii = 0;
    jj = 3 * segmentLen;
    while (jj < 4 * segmentLen) {
        array[ii] = pWaveletNega[jj];
        ii++;
        jj++;
    }
    temp4 = ECGMin(array, 4 * segmentLen);

    temp = (temp1 + temp2 + temp3 + temp4) / 4;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //% nega=-1*(nega<thnega*negaThreshold);%% 1/10
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float waveletNega[pLength] = {0};//zeros(1,pLength);
    ii = 0;
    while (ii < pLength) {
        if (pWaveletNega[ii] < temp * negaThreshold)
            waveletNega[ii] = -1;
        else {
            waveletNega[ii] = 0;
        }
        ii++;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //% interva=posi+nega;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float waveletP[pLength] = {0};//zeros(1,pLength);
    ii = 0;
    while (ii < pLength) {
        waveletP[ii] = waveletNega[ii] + waveletPosi[ii];
        ii++;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %  loca=find(interva);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int pWaveletLoc[20] = {0};//zeros(1,20);

    ii = 0;
    jj = 0;
    while (ii < pLength) {
        if (waveletP[ii] != 0) {
            pWaveletLoc[jj] = (int) ii;
            jj = jj + 1;
        }
        ii++;
    }

    //   %%%%%%%%%%%%%%%%%% test   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
//     loca=pWaveletLoc;
    //   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unsigned int pplen = jj;

    float negMax = maxInteger;
    float posiMax = minInteger;

    ii = 0;
    while (ii < pplen) {
        if (((pMj4[pWaveletLoc[ii]]) > posiMax) && (pMj4[pWaveletLoc[ii]]) > 0) {
            posiMax = pMj4[pWaveletLoc[ii]];
        } else {
            if (((pMj4[pWaveletLoc[ii]]) < negMax) && (pMj4[pWaveletLoc[ii]]) < 0) {
                negMax = pMj4[pWaveletLoc[ii]];
            }
        }
        ii++;
    }

    int locaLength = 0;
    //for i=1:length(loca)
    ii = 0;
    while (ii < pplen) {
        if (((pMj4[pWaveletLoc[ii]] / posiMax) > 0.096) || (pMj4[pWaveletLoc[ii]] / negMax) > 0.096) {
            locaLength = locaLength + 1;
        }

        if (((pMj4[pWaveletLoc[ii]]) > 0) && ((pMj4[pWaveletLoc[ii]] / posiMax) < 0.096)) {
            pWaveletLoc[ii] = 0;
        } else if (((pMj4[pWaveletLoc[ii]]) < 0) && ((pMj4[pWaveletLoc[ii]] / negMax) < 0.096)) {
            pWaveletLoc[ii] = 0;
        }

        if ((ii > 1) && (pWaveletLoc[ii - 1] > 100) && (pWaveletLoc[ii] != 0))
            if ((((pMj4[pWaveletLoc[ii - 1]]) * (pMj4[pWaveletLoc[ii]])) < 0) &&
                ((fabsf(pMj4[pWaveletLoc[ii - 1]] / pMj4[pWaveletLoc[ii]]) < 0.096))) {
                pWaveletLoc[ii - 1] = 0;
                locaLength = locaLength - 1;
            }
        ii++;
    }


    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int localResult[20] = {0};//zeros(1,locaLength);
    ii = 0;
    jj = 0;
    // for i=1:length(loca)
    while (ii < pplen) {
        if (pWaveletLoc[ii] > 0) {
            localResult[jj] = pWaveletLoc[ii];
            jj++;
        }
        ii++;
    }

    // loca=localResult;
    float diff[20] = {0};
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ii = 0;
    while (ii < 20) {
        if ((abs(localResult[ii] - localResult[ii + 1]) < modulusDist) && ((localResult[ii]) != 0) &&
            ((localResult[ii + 1]) != 0)) {
            diff[ii] = waveletP[localResult[ii]] - waveletP[localResult[ii + 1]];//%%%%%  17:02 20181204
        } else {
            diff[ii] = 0;
        }
        ii++;
    }

    //  % look for maximum pairs
    //  %loca2=find(diff==-2);
    //  %loca2=find((diff==-2)|(diff==2));
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 0;
    jj = 0;
    int waveletLocation[20] = {0};//zeros(1,1);
    while (ii < 20) {
        if ((diff[ii] == 2) || (diff[ii] == -2)) {
            waveletLocation[jj] = (int) ii;
            jj++;
        }
        ii++;
    }

    int pWaveSpecial = 0;
    if ((jj >= 1) && ((diff[0] == 2) || (diff[0] == -2))) {
        pWaveSpecial = 1;
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // % % negative maximum modulum pairs
    // % interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
    // % %positive maximum modulum pairs
    // % interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    int pResult[pLength] = {0};//zeros(1,dataLength);
    unsigned int pLen = jj;//length(waveletLocation);
    ii = 0;
    int m, k;
    while (ii < pLen) {
        jj = waveletLocation[ii];
        if ((jj != 0) || (pWaveSpecial == 1)) {
            k = localResult[jj];
            m = localResult[jj + 1];
            //% negative maximum modulum pairs
            pResult[k] = (int) waveletP[k];
            // % positive maximum modulum pairs
            pResult[m] = (int) waveletP[m];
        }
        ii = ii + 1;
    }

//   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (pLen >= 6)//%%&&(abs( pPeak)<(abs(ecgdata(Ppos))))
    {

        int Ppos;
        pWaveRoutine(Xpc, Ponset, Poffset, &Ppos);
        *pPeak = Xpc[Ppos];
        *pPeakPositionFirst = Ppos;
        if (Xpc[Ppos] > 0)
            *type = 1;
        else
            *type = 2;

        pSumCalculate(Xpc, *Ponset, *Poffset, Psum);
        *pPeakPositionSecond = 0;

        return;
    }

    int mark1 = 0;
    int mark2 = 0;
    int mark3 = 0;
    int mark2Bak = 0;

    i = onset - 129 + pError - 1;
    if (i < 1)
        i = 60;
    // end
    if (pLength < 2 && rr < BEAT_MS600)
        //   i=BEAT_MS60;
        //end
        if (rr < BEAT_MS400) {
            i = BEAT_MS140;
        }

    int PeakPosition = 0;
    float limit2 = 0.15f;
    int workFlag = 0;
    int Pleft;
    unsigned int begin;
    // float slope[100];
    int p[10] = {0};
    double radium;
    int pNumber = 0;



//%%%%%%%%%%%%%%%%%%%%% main route   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%% main route   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while ((i < (onset + pError * 2)) && (i < pLimit) && (pLen > 0))//~isempty(waveletLocation )))
    {
        if ((pResult[i] == -1) || (pResult[i] == 1)) {
            mark1 = (int) i;
            i = i + 1;
            while ((i < pLimit) && (pResult[i] == 0) && (mark2 == 0)) {
                i++;
            }

            if (mark2 == 0) {
                mark2 = (int) i;
            } else {
                temp = (float) mark1;
                mark1 = mark2;
                mark2 = (int) temp;
            }
            //  % modulus pair zero_crossing
            if ((pLen == 1) && (rr < BEAT_MS600)) {
                limit2 = 0.08f;
            }
            workFlag = 1;
            int Peak;
            int position;

            if ((((pMj4[mark1] * pMj4[mark2]) < 0) && (abs(mark2 - mark1) < BEAT_MS150)) &&
                ((fabsf(pMj4[mark2] / pMj4[mark1]) > limit2) &&
                 fabsf(pMj4[mark1] / pMj4[mark2]) > limit2))//%&&abs(onset-PeakPosition)<100)
            {
                mark3 = (int) roundf((fabsf(pMj4[mark2]) * (float) mark1 + (float) mark2 * fabsf(pMj4[mark1])) /
                                     (fabsf(pMj4[mark2]) + fabsf(pMj4[mark1])));
                Peak = mark3 - pError * 3;
                if (Peak < 0) {
                    Peak = 1;
                }
                k = Peak - 16;
                if (k < 0) {
                    k = 0;
                }
                temp = Xpc[k];
                position = k;

                while (k < (Peak + 10)) {
                    if ((pMj4[mark1] < 0) && pMj4[mark2] > 0)  // %negative maximum modulus pairs
                    {
                        if (Xpc[k] > temp) {
                            temp = Xpc[k];
                            position = k;
                        }
                    } else if ((pMj4[mark1] > 0) && pMj4[mark2] < 0)// %positive maximum modulus pairs
                    {
                        if (Xpc[k] < temp) {
                            temp = Xpc[k];
                            position = k;
                        }
                    }

                    k = k + 1;
                }
                // %%%%%%%%%%%
                PeakPosition = position;
                Pleft = mark1 - pError - 30;
                begin = Pleft - 5;

                if (((PeakPosition - BEAT_MS120) > 0) && (pLen == 1)) {
                    begin = PeakPosition - BEAT_MS120;
                }
                if (begin < 1) {
                    begin = 6;
                }
                if ((PeakPosition - BEAT_MS16) < 0) {
                    continue;
                }
                *Ponset = localMaxAlgorithm(Xpc, (int) begin, PeakPosition - BEAT_MS16);
                // for many steps before %%%%%%
                // /*
                if ((abs(PeakPosition - *Ponset) < 20) &&
                    abs((int) onset - PeakPosition) < 106)//<200ms  16:32 11/08 2018  ypz
                {
                    m = 0;
                    float maxSlopeI = 0;
                    float minSlopeI = 0;
                    k = (int) begin;

                    float maxSlope = Xpc[k + 1] - Xpc[k];
                    float minSlope = Xpc[k + 1] - Xpc[k];
                    float currentSlope = maxSlope;
                    float lastSlope = maxSlope;
                    while (k < (PeakPosition - 1)) {
                        //slope[k]=Xpc[k+1]-Xpc[k] ;
                        lastSlope = currentSlope;
                        currentSlope = Xpc[k + 1] - Xpc[k];
                        // if(slope[k] > maxSlope)
                        if (currentSlope > maxSlope) {
                            maxSlope = currentSlope;//slope[k] ;    //  max slope
                            maxSlopeI = (float) k;
                        } else if (currentSlope < minSlope) {
                            minSlope = currentSlope;//slope[k];  // min slope
                            minSlopeI = (float) k;
                        }

                        // if((k-1)<1)
                        //   if(((k-1)>=1)&&(slope[k]*slope[k-1])<0)
                        if (((k - 1) > 0) && (currentSlope * lastSlope) < 0) {
                            p[m] = k;
                            m = m + 1;                   // zero crosing
                        }
                        k = k + 1;
                    }

                    int found = 0;
                    m = m - 1;
                    while ((m > 0) && (found > 0))  // meet the  requirements ,nearset point
                    {
                        if ((abs(PeakPosition - p[m]) >= 20) && (p[m] != 0) && (k < *Ponset)) {
                            *Ponset = p[m];
                            found = 1;
                        }
                        m = m - 1;
                    }
                }// end 	if((fabs(PeakPosition-*Ponset)<20)&&fabs(onset-PeakPosition)<106)//<
// */
                //   %%% for next search
                /*
                    if((found==0)&&(m!=0))   //      %% monotonic rise/down  line
                           {
                             if((pMj4[mark1]<0)&&(pMj4[mark2]>0))
                               {
                                 k=*Ponset;
                                 int pFound=0;
                                 while((k>=Pleft)&&(pFound==0))
                                  {
                                     if((slope[k]<maxSlope/20)&&(fabs(PeakPosition-k)>=20))
                                      {
                                       *Ponset=k;
                                        pFound=1;
                                      }
                                    k=k-1;
                                  }
                                   if(pFound==0)
                                   {
                                   *Ponset=k;  //final result 50ms
                                   }
                               }

                           }// end  if((found==0)&&(m!=0))
                        }//end   if((abs(PeakPosition-*Ponset)<20)&&fabs(onset-PeakPosition)<106)

*/

                if (PeakPosition < 0) {
                    continue;
                }
                // %%%%%%%%%%%%%%%%%%%%%%
                int pRight = mark2;
                if ((pRight > (pLimit)))
                    pRight = pLimit - 10;
                begin = (unsigned int) pRight;
                if (begin > onset)
                    begin = onset;

                *Poffset = localMaxAlgorithm(Xpc, PeakPosition + BEAT_MS40, (int) begin);
                if ((*Poffset < (mark2 - pError)) && ((mark2 - pError + 25) <
                                                      (int) onset))//%%&&(ecgdata(PeakPosition)>(ecgdata(mark2-pError)+ecgdata(mark1-pError))/2)
                {

                    if ((*Poffset < 0) || (mark2 - pError + 25) < 0) {
                        continue;
                    }
                    *Poffset = localMaxAlgorithm(Xpc, *Poffset, mark2 - pError + 25);
                }

                if ((abs(*Poffset - *Ponset) >= 24) && ((onset - *Ponset) < 160) &&
                    (abs(PeakPosition - *Ponset) >= 13)) {
                    *Ponset = floor(*Ponset);
                    *Poffset = floor(*Poffset);
                    if (*Ponset < 0)
                        *Ponset = 1;

                    if (*Poffset < 0)
                        *Poffset = 1;

                    double distance = sqrt((double) (*Poffset - *Ponset) * (double) (*Poffset - *Ponset) +
                                           ((double) Xpc[*Poffset] - (double) Xpc[*Ponset]) *
                                           ((double) Xpc[*Poffset] - (double) Xpc[*Ponset]));
                    int x0 = PeakPosition;
                    int x1 = *Ponset;
                    int x2 = *Poffset;
                    float y0 = Xpc[x0];
                    float y1 = Xpc[x1];
                    float y2 = Xpc[x2];
                    float A = y2 - y1;
                    float B = (float) (x1 - x2);
                    float C = (float) x2 * y1 - (float) x1 * y2;
                    double d = fabs(((double) A * (double) x0 + (double) B * (double) y0 + (double) C) /
                                    sqrt((double) A * (double) A + (double) B * (double) B));
                    radium = d / distance;
                }

                if ((radium < limit)) {
                    workFlag = 0; // default
                } else {          //Fs<173e-3  193e-3
                    if (((double) abs(*Poffset - *Ponset) / Fs < 220e-3) && (abs(PeakPosition - *Ponset) >= 13) &&
                        abs(PeakPosition - *Poffset) >= 11)//%%40ms
                    {

                        Presult[pNumber].onset = *Ponset;
                        Presult[pNumber].Peak = PeakPosition;
                        Presult[pNumber].offset = *Poffset;
                        Presult[pNumber].weight = Xpc[PeakPosition];
                        Presult[pNumber].coeff = (float) radium;
                        mark2Bak = mark2;

                        if (Xpc[PeakPosition] < (Xpc[*Ponset] + Xpc[*Poffset]) / 2) {
                            Presult[pNumber].type = 2;//type_2;
                        } else {
                            if (Xpc[PeakPosition] > (Xpc[*Ponset] + Xpc[*Poffset]) / 2) {
                                Presult[pNumber].type = 1;//type_1;
                            }
                        }
                        pNumber = pNumber + 1;
                    }
                }// end if((radium<limit)||abs(Mj4(mark2)/Mj4(mark1))<limit||abs(Mj4(mark1)/Mj4(mark2))<limit)
            } //end   if((((pMj4[mark1]*pMj4[mark2])<0)&&(fabs(mark2-mark1)<BEAT_MS150))&&((fabs(pMj4[mark2]/pMj4[mark1])>limit2)&&f
        } //end   if((pResult[i]==-1)||(pResult[i]==1))

        i = i + 1;
    }//end    while((i<(onset+pError*2))&&(i<pLimit)&&(pLen>0))//~isempty(waveletLocation )))

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ((mark2Bak != 0))//%%&&(mark1Bak~=0))
        mark2 = mark2Bak;
    //end

    if ((workFlag == 1) && (pNumber != 0)) {
        i = mark2 + 1;
        int Flag = 0;
        while ((i < (onset)) && (Flag == 0)) {
            if (pMj4[i] != 0) {
                if ((pMj4[i] > 0) && (pMj4[i] / pMj4[mark2] > 1)) {
                    temp = (float) i;
                    Flag = 1;
                }
            }
            i = i + 1;
        }
        if (Flag == 1) {
            int pEnd = (int) (temp + 5);
            if (pEnd > pLimit)
                pEnd = pLimit;

            if (pNumber > 1) {
                result.peak2 = Presult[pNumber - 1].offset;
                if (((temp - 20) > 0) && (pEnd) > 0) {
                    Presult[pNumber - 1].offset = localMaxAlgorithm(Xpc, (int) (temp - 20), pEnd);//PmaxPosition;
                }

            }
        }
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%thirdWeight

    if ((Presult[0].onset) >= BEAT_MS20) {
        float sum = 0;
        k = 0;
        while (k < BEAT_MS20) {
            sum = sum + (Xpc[Presult[0].onset - k]);
            k = k + 1;
        }
        temp = sum / BEAT_MS20;
        if (fabsf(temp) < fabsf(isoelectricLevel))
            isoelectricLevel = temp;

        if ((fabsf(isoelectricLevel / beat[FIDMARK]) > 1.0f / 20) && (beat[FIDMARK] > 0))
            isoelectricLevel = 0;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch (pNumber) {
        case 0:
            if ((PeakPosition != 0)) {
                if (((onset - *Ponset) > 110) && ((PeakPosition - *Ponset) < 20) && abs(Xpc[PeakPosition] > 0)) {
                    result.type = 5;//type_5;
                } else {
                    int start = FIDMARK - 0.25 * Fs;
                    int finish = FIDMARK - 0.1 * Fs;
                    float temp[76];
                    bufferCopy(Xpc, start, finish, temp);

                    float returnMax = 0;
                    int position = 0;
                    ECGMaxPosition(Xpc, (int) (0.15 * Fs), &returnMax, &position);

                    //%%%%%%%%%%%%%%
                    if (*Ponset > 0 && fabsf(returnMax) > 3 && *Poffset > 0 &&
                        abs((int) onset - *Poffset) < BEAT_MS120) {
                        if (returnMax > 5) {
                            result.type = 1;//type_1;
                            result.onset = 30;
                            result.peak1 = position + start;
                            result.offset = 120;
                        } else {   // junctional inversive P wave
                            result.type = 2;//type_2;
                            result.onset = *Ponset;
                            result.peak1 = position + start;
                            result.offset = *Poffset;
                        }
                    }
                }
            }

            break;

//  ###############################################
        case 1:
            result.onset = Presult[0].onset;
            result.peak1 = Presult[0].Peak;
            result.peak2 = 0;
            result.offset = Presult[0].offset;
            result.type = Presult[0].type;
            result.weight = Presult[0].weight;
            if (((Presult[0].offset - Presult[0].onset) <= BEAT_MS110) && ((onset - Presult[0].onset) > BEAT_MS210)) {
                result.type = 6;//type_6;
            } else {
                if ((result.peak1 != 0) && (result.onset != 0)) {
                    m = (int) fabsf((Xpc[result.peak1] - Xpc[result.onset]) / (float) (result.peak1 - result.onset));
                    if (((m > 1.5) && ((onset - result.onset) > 110)) ||
                        (((onset - *Ponset) > BEAT_MS220) && ((onset - PeakPosition)) > BEAT_MS200)) {
                        result.type = type_5;
                    }
                }
                if (((onset - result.onset) > 110) && ((onset - result.peak1) > 100)) {
                    result.type = type_5;
                }
            }

            break;

            //  ###############################################
        case 2: {
            // %% mathematic morphology compare
            float ration;
            if (((Presult[0].coeff)) > ((Presult[1].coeff)))
                ration = ((Presult[1].coeff)) / (((Presult[0].coeff)));
            else
                ration = ((Presult[0].coeff)) / (((Presult[1].coeff)));

            float firstWeight = Presult[0].weight - isoelectricLevel;
            float secondWeight = Presult[1].weight - isoelectricLevel;
            float heightRatio;

            if (fabsf(firstWeight) > fabsf(secondWeight)) {
                heightRatio = fabsf(firstWeight / secondWeight);
            } else
                heightRatio = fabsf(secondWeight / firstWeight);


            if ((firstWeight > fabsf(isoelectricLevel * 2)) && (secondWeight > fabsf(isoelectricLevel))) {
                if (((Presult[1].type) == type_1) && (Presult[0].type) == type_1) {
                    result.onset = Presult[0].onset;
                    result.peak1 = Presult[0].Peak;
                    result.peak2 = Presult[1].Peak;
                    result.offset = Presult[1].offset;
                    result.type = type_1;
                } else {
                    result.onset = Presult[0].onset;
                    if ((firstWeight) > (secondWeight)) {
                        result.peak1 = Presult[0].Peak;
                    } else {
                        result.peak1 = Presult[1].Peak;
                    }
                    result.offset = Presult[1].offset;
                    result.type = type_1;

                    // %%%%%%%%%%% Ectopic beat %%%%%%%%%%%%%%%%%%%
                    if ((onset - Presult[0].onset) > 120)
                        result.type = type_5;
                }
            } else {
                if (((Presult[1].type) == type_1) && (Presult[0].type) == type_2) {
                    if ((abs((int) onset - Presult[1].onset) < BEAT_MS120) &&
                        (fabsf(secondWeight / (firstWeight)) > 1.0)) {
                        result.onset = Presult[1].onset;
                        result.peak1 = Presult[1].Peak;
                        result.offset = Presult[1].offset;
                        result.type = Presult[1].type;
                    } else {
                        if (ration > 0.7 &&
                            (float) abs(Presult[1].offset - Presult[0].onset) / Fs < 170e-3)// %% mix one peak
                        {
                            result.onset = Presult[0].onset;
                            result.peak1 = Presult[0].Peak;
                            result.peak2 = Presult[1].Peak;
                            result.offset = Presult[1].offset;
                            result.type = type_4; //  double wave  type 3
                        } else if (fabsf(Presult[0].coeff) > fabsf(Presult[1].coeff)) {
                            {
                                result.onset = Presult[0].onset;
                                result.peak1 = Presult[0].Peak;
                                result.offset = Presult[0].offset;
                                result.type = Presult[0].type;
                            }
                            if (((Presult[0].offset - Presult[0].onset) <= BEAT_MS110) &&
                                ((onset - Presult[0].onset) > BEAT_MS210)) {
                                result.type = type_6;   //  I AVB P wave
                            } else {
                                if (((onset - result.onset) > 110) && ((onset - result.peak1) > 100)) {
                                    result.type = type_5; // ectopic P wave

                                }
                            }
                        } else {
                            // % [ Angle ] =findDegree( ecgdata,Presult(2).onset,Presult(2).Peak,Presult(2).offset )firstWeight
                            if (fabsf(firstWeight / secondWeight) > 1.2) {
                                result.onset = Presult[0].onset;
                                result.peak1 = Presult[0].Peak;
                                result.offset = Presult[0].offset;
                                result.type = Presult[0].type;
                            } else {
                                result.onset = Presult[1].onset;
                                result.peak1 = Presult[1].Peak;
                                result.offset = Presult[1].offset;
                                result.type = Presult[1].type;
                            }
                        }
                    }

                } else  //   *********
                {
                    if (((Presult[1].type) == 1) && (Presult[0].type) == 1)  //  mode 1,1
                    {
                        if (fabsf(Presult[1].weight) > fabsf(Presult[0].weight)) {
                            result.onset = Presult[1].onset;
                            result.peak1 = Presult[1].Peak;
                            result.offset = Presult[1].offset;
                            result.type = Presult[1].type;
                        } else {
                            {
                                result.onset = Presult[0].onset;
                                result.peak1 = Presult[0].Peak;
                                result.offset = Presult[0].offset;
                                result.type = Presult[0].type;
                            }

                            // %  if(((onset-result.onset)>110)&&((result.peak1-result.onset)<30)&&abs(ecgdata(result.peak1)>0))
                            if (((onset - result.onset) > 110) && ((onset - result.peak1) > 100)) {
                                result.type = type_5;
                            }

                        }
                    } else {
                        if (((Presult[1].type) == type_2) && (Presult[0].type) == type_2) // mode 2,2
                        {
                            if ((fabsf(Presult[1].coeff)) < (fabsf(Presult[0].coeff)) ||
                                (abs((int) onset - Presult[1].onset) < 109)) {
                                result.onset = Presult[1].onset;
                                result.peak1 = Presult[1].Peak;
                                result.offset = Presult[1].offset;
                                result.type = Presult[1].type;
                            } else {
                                result.onset = Presult[0].onset;
                                result.peak1 = Presult[0].Peak;
                                result.offset = Presult[0].offset;
                                result.type = Presult[0].type;
                            }
                        }
                    }
                }

            }

            break;
        }

            //  ###############################################

        case 3: {
            float firstWeight = Presult[0].weight - isoelectricLevel;
            float secondWeight = Presult[1].weight - isoelectricLevel;
            float thirdWeight = Presult[2].weight - isoelectricLevel;


            int rationCoeff = fabsf(thirdWeight / secondWeight) > 2 && fabsf(thirdWeight / firstWeight) > 2;
            if (rationCoeff == 1) {
                *Ponset = Presult[2].onset;
                *pPeak = Presult[2].weight;
                *Poffset = Presult[2].offset;
                *pPeakPositionFirst = Presult[2].Peak;
                if (*pPeak > 0)
                    *type = type_1;
                else
                    *type = type_2;

                if ((*Ponset > 0) && (*Poffset > 0)) {
                    pSumCalculate(Xpc, *Ponset, *Poffset, Psum);
                    pPeakPositionSecond = 0;
                    return;
                }
            }

            if (fabsf(firstWeight / secondWeight) > 2 && fabsf(firstWeight / thirdWeight) > 2) {
                result.onset = Presult[0].onset;
                result.peak1 = Presult[0].Peak;
                result.offset = Presult[0].offset;
                result.type = Presult[0].type;

            } else {
                if (abs((int) onset - Presult[1].onset) < 118) {
                    float coef = fabsf(thirdWeight / secondWeight); //  201811129nmark
                    int tempTtpe =
                            (Presult[0].type == type_2) && (Presult[1].type == type_1) && (Presult[2].type == type_2);
                    if ((abs((int) onset - Presult[2].onset) < BEAT_MS120) && (Presult[2].type == type_2) &&
                        (coef >= 1) &&
                        (fabsf(thirdWeight / beat[FIDMARK]) > 0.02) && (tempTtpe == 0)) {
                        result.onset = Presult[2].onset;
                        result.peak1 = Presult[2].Peak;
                        result.offset = Presult[2].offset;
                        result.type = Presult[2].type;
                    } else {
                        if ((coef > 0.5) && (coef < 2)) {
                            result.onset = Presult[1].onset;
                            result.peak1 = Presult[1].Peak;
                            result.peak2 = Presult[2].Peak;
                            result.offset = Presult[2].offset;
                            result.type = type_4;
                        } else if (fabsf(thirdWeight) > fabsf(secondWeight)) {
                            result.onset = Presult[2].onset;
                            result.peak1 = Presult[2].Peak;
                            result.offset = Presult[2].offset;
                            result.type = Presult[2].type;
                        } else {
                            result.onset = Presult[1].onset;
                            result.peak1 = Presult[1].Peak;
                            result.offset = Presult[1].offset;
                            result.type = Presult[1].type;
                        }

                    }// end    if(coef>0.5)&&(coef<2)
                } else {
                    temp = 0;
                    int number = 0;
                    if (fabsf(Presult[0].weight - isoelectricLevel) > fabsf(Presult[1].weight - isoelectricLevel)) {
                        temp = fabsf(Presult[0].weight - isoelectricLevel);
                        number = 1;
                    }
                    if (fabsf(Presult[1].weight - isoelectricLevel) > temp) {
                        temp = fabsf(Presult[1].weight - isoelectricLevel);
                        number = 2;
                    }
                    if (fabsf(Presult[2].weight - isoelectricLevel) > temp) {
                        temp = fabsf(Presult[2].weight - isoelectricLevel);
                        number = 3;
                    }


                    switch (number) {
                        case 1: {
                            float pTvalue = fabsf((Presult[1].weight - isoelectricLevel) /
                                                  (Presult[number].weight - isoelectricLevel));
                            if ((pTvalue > 0.6) && (abs((int) onset - Presult[number].onset) < 105)) {
                                result.onset = Presult[number].onset;
                                result.peak1 = Presult[number].Peak;
                                result.peak2 = 0;
                                result.offset = Presult[number].offset;
                                result.type = type_1;
                                if (((onset - result.onset) > 110) && ((onset - result.peak1) > 100)) {
                                    result.type = type_5;
                                }
                            } else {
                                if (((Presult[1].type) == type_1) && (Presult[1].type == type_2))
                                    if (((onset - Presult[number].onset) > 110) &&
                                        ((onset - Presult[number].Peak) > 100)) {
                                        result.onset = Presult[number].onset;
                                        result.peak1 = Presult[number].Peak;
                                        result.offset = Presult[number].offset;
                                        result.type = type_5;
                                    } else {
                                        result.onset = Presult[0].onset;
                                        result.peak1 = Presult[0].Peak;
                                        result.peak2 = Presult[1].Peak;
                                        result.offset = Presult[1].offset;
                                        result.type = type_3;
                                    }
                                else {
                                    result.onset = Presult[0].onset;
                                    result.peak1 = Presult[0].Peak;
                                    result.peak2 = Presult[1].Peak;
                                    result.offset = Presult[1].offset;
                                    result.type = type_4;
                                }
                            }

                        }
                            break;


                        case 2: {

                            float pTemp = fabsf(
                                    (Presult[1].weight - isoelectricLevel) / (Presult[0].weight - isoelectricLevel));
                            float pTemp2 = fabsf(
                                    (Presult[1].weight - isoelectricLevel) / (Presult[2].weight - isoelectricLevel));
                            float pTemp4 = ((Presult[1].coeff) / (Presult[2].coeff));
                            int pType = ((Presult[0].type == type_2) && (Presult[1].type == type_1) &&
                                         (Presult[2].type == type_2) &&
                                         abs((int) onset - Presult[1].onset) < 105);// pr interval within about 200ms
                            if ((pType == type_1) && (pTemp > 1.45) && (pTemp2 >= 1.85)) {
                                result.onset = Presult[1].onset;
                                result.peak1 = Presult[1].Peak;
                                result.offset = Presult[1].offset;
                                result.type = type_1;
                            } else {
                                if ((pTemp2 >= 1) && (pTemp2 <= 2)) {
                                    if (((Presult[2].Peak - Presult[1].offset) < 5) && (pTemp4 > 0.6))//  adjcent
                                    {
                                        if ((Presult[2].type) != Presult[1].type) {
                                            result.onset = Presult[1].onset;
                                            result.peak1 = Presult[1].Peak;
                                            result.peak2 = Presult[2].Peak;
                                            result.offset = Presult[2].offset;
                                            if (Presult[1].type == 1) {
                                                result.type = type_3;
                                            } else {
                                                result.type = type_4;
                                            }
                                        } else {//  double peaks
                                            result.onset = Presult[1].onset;
                                            result.peak1 = Presult[1].Peak;
                                            result.peak2 = Presult[2].Peak;
                                            result.offset = Presult[2].offset;
                                            result.type = Presult[1].type;
                                        }
                                    } else
                                        //    not adjcent
                                    {
                                        if ((Presult[2].coeff) > (Presult[1].coeff)) {
                                            result.onset = Presult[2].onset;
                                            result.peak1 = Presult[2].Peak;
                                            result.offset = Presult[2].offset;
                                            result.type = Presult[2].type;
                                        } else {
                                            result.onset = Presult[1].onset;
                                            result.peak1 = Presult[1].Peak;
                                            result.offset = Presult[1].offset;
                                            result.type = Presult[1].type;
                                        }
                                    }

                                } else {
                                    result.type = type_5;

                                }

                            }

                        }
                            break;

                        case 3: {

                            if (((fabsf(thirdWeight / secondWeight)) > 1) &&
                                ((fabsf(thirdWeight / secondWeight)) < 1.5)) {
                                if (((Presult[1].offset - Presult[1].onset) <= BEAT_MS110) &&
                                    ((onset - Presult[1].onset) > BEAT_MS210)) {
                                    result.onset = Presult[1].onset;
                                    result.peak1 = Presult[1].Peak;
                                    result.offset = Presult[1].offset;
                                    result.type = type_6; // might be Ectopic beat
                                } else {
                                    result.onset = Presult[2].onset;
                                    result.peak1 = Presult[2].Peak;
                                    result.offset = Presult[2].offset;
                                    result.type = Presult[2].type;
                                }
                            } else if (Presult[1].offset - Presult[1].onset <= BEAT_MS110
                                       && onset - Presult[1].onset > BEAT_MS210) {
                                result.onset = Presult[1].onset;
                                result.peak1 = Presult[1].Peak;
                                result.offset = Presult[1].offset;
                                result.type = type_6;
                            } else if (fabsf(Presult[1].weight / Presult[number].weight) < 0.6) {
                                result.onset = Presult[2].onset;
                                result.peak1 = Presult[2].Peak;
                                result.peak2 = 0;
                                result.offset = Presult[2].offset;
                                result.type = type_1;

                                if (((onset - result.onset) > 110) && ((onset - result.peak1) > 100)) {
                                    result.type = type_5;
                                }
                            } else {
                                if (((Presult[1].type) == type_1) && (Presult[2].type == type_2)) {
                                    result.onset = Presult[1].onset;
                                    result.peak1 = Presult[1].Peak;
                                    result.peak2 = Presult[2].Peak;
                                    result.offset = Presult[2].offset;
                                    result.type = type_3;
                                } else {
                                    result.onset = Presult[1].onset;
                                    result.peak1 = Presult[1].Peak;
                                    result.peak2 = Presult[2].Peak;
                                    result.offset = Presult[2].offset;
                                    result.type = type_4;
                                }
                            }
                            // end  (((onset-Presult(1).onset)>110)&&((onset-Presult(1).Peak)>100))
                            //  end  //(((Presult(2).offset-Presult(2).onset)<=BEAT_MS110)&&((onset-Presult(2).onset)>BEAT_MS210))
                            // end

                        }
                        default:;
                            break;
                    }

                }

            }//end if ((abs(firstWeight/secondWeight))>2 )&&((abs(firstWeight/thirdWeight))>2)

            break;
        }

//  ###########  case 4:  ####################################
        case 4: {
            if ((Presult[2].type == type_2) && (Presult[3].type == type_1)) {
                *Ponset = Presult[2].Peak;
                int Ppos = Presult[3].Peak;
                *Poffset = Presult[3].offset;
                *type = type_1;
                *pPeakPositionFirst = Ppos;
                pSumCalculate(Xpc, *Ponset, *Poffset, Psum);
                pPeakPositionSecond = 0;
                return;
            }

            float pWeightArray[4] = {0};
            float pIntervalArray[4] = {0};
            pWeightArray[0] = Presult[0].weight - isoelectricLevel;
            pWeightArray[1] = Presult[1].weight - isoelectricLevel;
            pWeightArray[2] = Presult[2].weight - isoelectricLevel;
            pWeightArray[3] = Presult[3].weight - isoelectricLevel;

            pIntervalArray[0] = (float) onset - (float) Presult[0].onset;
            pIntervalArray[1] = (float) onset - (float) Presult[1].onset;
            pIntervalArray[2] = (float) onset - (float) Presult[2].onset;
            pIntervalArray[3] = (float) onset - (float) Presult[3].onset;

            float pWeightMax = pWeightArray[0];
            int pWaveInterval = 1;
            k = 0;
            while (k < 4) {
                if ((pIntervalArray[k] < BEAT_MS210) && (pWeightArray[k] > pWeightMax)) {
                    pWeightMax = pWeightArray[k];
                    pWaveInterval = k;
                }
                k = k + 1;
            }

            result.onset = Presult[pWaveInterval].onset;
            result.peak1 = Presult[pWaveInterval].Peak;
            result.offset = Presult[pWaveInterval].offset;
            result.type = Presult[pWaveInterval].type;

            if ((abs((int) onset - Presult[3].onset) < BEAT_MS120) && (Presult[3].type == type_2) &&
                ((fabsf(pWeightArray[3] / pWeightMax)) > 0.5)) {
                result.onset = Presult[3].onset;
                result.peak1 = Presult[3].Peak;
                result.offset = Presult[3].offset;
                result.type = Presult[3].type;
            }

        }
            break;

//  ###########  case 5:  ####################################

        case 5: {

            float pWeight2 = Presult[1].weight - isoelectricLevel;
            float pWeight3 = Presult[2].weight - isoelectricLevel;
            float pWeight4 = Presult[3].weight - isoelectricLevel;

            if ((pWeight2 > 0) && (pWeight3 > 0) && (pWeight4 > 0)) {
                result.onset = Presult[1].onset;
                int position = threeDataMax(pWeight2, pWeight3, pWeight4);
                result.peak1 = Presult[position].Peak;
                result.offset = Presult[3].offset;
                result.type = type_1;
            }
        }
            break;

//  ###########  default ####################################
        default: {
            int pPonset;
            int pPoffset;
            int Ppos;
            pWaveRoutine(Xpc, &pPonset, &pPoffset, &Ppos);
            *Ponset = pPonset;
            *pPeak = Xpc[Ppos];
            *Poffset = pPoffset;
            *pPeakPositionFirst = Ppos;
            if (Xpc[Ppos] > 0) {
                *type = type_1;
            } else {
                *type = type_2;
            }
            pSumCalculate(Xpc, *Ponset, *Poffset, Psum);
            pPeakPositionSecond = 0;
            return;
        }

    }// end   switch(pNumber)

//  ###########  end switch  ####################################

    if ((result.peak1) != 0) {
        *pPeakPositionFirst = result.peak1;
    } else
        *pPeakPositionFirst = 0;

    if ((result.peak2) != 0) {
        *pPeakPositionSecond = result.peak2;
    } else
        *pPeakPositionSecond = 0;

    *Ponset = result.onset;


    if ((result.peak1 != 0) && (result.peak2 != 0))
        if (fabsf(Xpc[result.peak2]) > fabsf(Xpc[result.peak1])) {
            *pPeak = beat[result.peak2];
            *pPeakPositionFirst = result.peak2;

        } else {
            *pPeak = beat[result.peak1];
            *pPeakPositionSecond = result.peak1;

        }

    else if (result.peak1 != 0) {
        *pPeak = beat[result.peak1];
        *pPeakPositionFirst = result.peak1;
    }

    *Poffset = result.offset;
    *type = result.type;
    pSumCalculate(Xpc, *Ponset, *Poffset, Psum);

    if ((*pPeak == 0) && (result.type == type_5) && (PeakPosition != 0))
        *pPeak = beat[PeakPosition];

    if (*pPeakPositionFirst == 156)
        ii = 0;

}
