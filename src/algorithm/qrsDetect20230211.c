
#include "qrsDetect20230211.h"
#include <math.h>
#include "common_util.h"
#define   BEAT_MS50    25
#define   BEAT_MS280   140
#define   BEAT_MS110   55
#define   FIDMARK      200


#define EPSILON 0.000001

//qrsRoute::qrsRoute()
//{
////   int kk=0;
////   for(kk=0;kk<300;kk++)
////   {
////     result2[kk]=0;

////   }


//}


//qrsRoute::~qrsRoute()
//{




//}


//  float qrsRoute::qrsHFNoiseCheck(float*beat,float qrsMin,float qrsMax)
float qrsHFNoiseCheck(float *beat, float qrsMin, float qrsMax) {


    int AVELENGTH = BEAT_MS50;

    float returnValue = 0;
    float maxNoiseAve = 0;
    float sum = 0;
    int avePtr = 0;
    float aveBuff[BEAT_MS50] = {0};


    int i = FIDMARK - BEAT_MS280;
    while (i < FIDMARK + BEAT_MS280) {
        sum = sum - aveBuff[avePtr];
        aveBuff[avePtr] = fabsf(beat[i] - (beat[i - BEAT_MS10] * 2) + beat[i - 2 * BEAT_MS10]);
        sum = sum + aveBuff[avePtr];
        avePtr = avePtr + 1;
        if (avePtr == AVELENGTH)
            avePtr = 0;
        //
        if ((i < (FIDMARK - BEAT_MS50)) || (i > (FIDMARK + BEAT_MS110)))
            if (sum > maxNoiseAve)
                maxNoiseAve = sum;
        // end
        // end

        i = i + 1;
    }

    if ((qrsMax - qrsMin) >= 4) {
        returnValue = ((maxNoiseAve * (50 / (float) AVELENGTH)) / ((qrsMax - qrsMin) / 4));
        return returnValue;
    } else {
        returnValue = 0;
        return returnValue;
    }

}



// /* initialize the static varies and array*/
//int  qrsRoute::r=0;
////float qrsRoute::TH =0.2125;//%0.3125;%0.2125;
////%   nodeInformation(r).qrsWidth=0;
////qrsFeatures(r).rr=0;
////qrsFeatures(r).qrsWidth=0;
//float qrsRoute:: rrArray[5]={0};
//float qrsRoute:: px[8]={0};
//float qrsRoute::py[8]={0};
////px=zeros(1,8);%[0.0 0.0 0.0 0.0 0.0 0.0 0.0];
////py=zeros(1,8);%[0.0 0.0 0.0 0.0 0.0 0.0 0.0];
//int qrsRoute::lastPeakIndex=0;
//int qrsRoute::adjcentPeak=0;
//int qrsRoute::lastPosition=0;
//int qrsRoute:: currentPosition1=0;
//float qrsRoute::rrbuf[8]={0};
//float qrsRoute::noise[8]={0};
//int qrsRoute::qrscount[maxQueueLength]={0};
////int qrsRoute::rrBuffer[maxQueueLength]={0};

//int qrsRoute::buffer[maxQueueLength]={0};
//SearchBackRRcount=0;
////ModNumber=maxQueueLength;
//int qrsRoute::RRcount=0;
//float qrsRoute::max = 0;
//int qrsRoute::timeSinceMax = 0;
//float qrsRoute::lastDatum=0 ;
////MIN_PEAK_AMP=10;
//float qrsRoute::PeakPosition[maxQueueLength] ={0};//zeros(1,maxQueueLength);
//float qrsRoute::lastPeak= MIN_PEAK_AMP;
////tempBeat=zeros(1,500);
//float qrsRoute::lostPeakArray[6]={0};//zeros(1,6);
//int qrsRoute::lostPeakPointer=0;
//float qrsRoute::lostPeakTemp=0;
//int qrsRoute::lostPeakOutputPointer=0;

//float qrsRoute::monintorArray[6]={0};//zeros(1,6);
//int qrsRoute::monitoringPointer=0;
////%%%%%%%%%%%%%%%%%%% moving window integrator Data Table %%%%%

////for MovingWindowSumPtr=1:1:WINDOW_WIDTH%%130%%%160
//float qrsRoute::data[WINDOW_WIDTH]={0};
////

//float qrsRoute::firstHeight=0;
//float qrsRoute::secondHeight=0;
////%  secondHeightPosition=0;
//float qrsRoute::MovingWindowSum=0;
//int qrsRoute::MovingWindowSumPtr=1;
////for i=1:1:10%60
//float qrsRoute::x_derv[10]={0};
////

////for i = 1:1:8
////%%noise(i) = 0 ;  %/* Initialize noise buffer */
////rrbuf(i)= MS1000 ;%/* and R-to-R interval buffer. */
//float qrsRoute::qrsbuf[8]={0};
//float qrsRoute::rsetBuff[8]={MS1000,MS1000,MS1000,MS1000,MS1000,MS1000,MS1000,MS1000};
////
//float qrsRoute::det_thresh=0;
//int qrsRoute::rsetCount = 1 ;
//float qrsRoute::nmean=0;
//float qrsRoute::qmean=0;
//float qrsRoute::rrmean =0;
//int qrsRoute::sbloc=0;
//int qrsRoute::QrsDelay = 0 ;
//float qrsRoute::lastNewPeak=0;
//float qrsRoute::currentNewPeak=0;

//int qrsRoute::smallCount=0;
//int qrsRoute::qpkcnt = 0;  //%%%%%%%%%%%%%%%%%%%
//int qrsRoute::count = 0;
//float qrsRoute::sbpeak = 0 ;
//int qrsRoute::initBlank = 0;
//float qrsRoute::initMax =0;
//int qrsRoute::preBlankCnt = 0;
//int qrsRoute::index=0;
//int qrsRoute::RRInterval=0;
//float qrsRoute::tempPeak=0 ;
//float qrsRoute::newPeak=0;
//float qrsRoute:: ECGBuffer[maxQueueLength]={0};
//int qrsRoute::lastNewPeakPosition=0;
//int qrsRoute::currentNewPeakPosition=0;
//int qrsRoute::sbcount=MS1500 ;
////float qrsRoute::result2[500]={0};


//static float PeakPosition[maxQueueLength];
//static int   qrscount[maxQueueLength];
//static float ECGBuffer[maxQueueLength];
////int qrsRoute::qrsDetect(float inputdata)//,int peakNumber,float*ECGBuffer)


// **************************************************************************************
// *******************  Main Route  *****************************************************
// **************************************************************************************

//int qrsDetect2023(float inputdata,int multiLeadNumber,float *ECGBuffer,int init)//
int qrsDetectMultiLeads(qrsDetectMultiLeads_CTX *ctx,float inputdata, int multiLeadNumber, int init) {

/* initialize the static varies and array*/
//  *****************************************************
    int *indexMulti = ctx->indexMulti; // static int indexMulti[leadsNumber];
    float *tempPeakMulti = ctx->tempPeakMulti; // static float tempPeakMulti[leadsNumber];
    float (*x_dervMulti)[10] = ctx->x_dervMulti; //static float x_dervMulti[leadsNumber][10];
    float *MovingWindowSumMulti = ctx->MovingWindowSumMulti; // static float MovingWindowSumMulti[leadsNumber];
    int *MovingWindowSumPtrMulti = ctx->MovingWindowSumPtrMulti; // static int MovingWindowSumPtrMulti[leadsNumber];
    float (*dataMulti)[40] = ctx->dataMulti; //static float dataMulti[leadsNumber][40];
    int *timeSinceMaxMulti = ctx->timeSinceMaxMulti; // static int timeSinceMaxMulti[leadsNumber];
    float *lastDatumMulti = ctx->lastDatumMulti; // static float lastDatumMulti[leadsNumber];
    float *maxMulti = ctx->maxMulti; // static float maxMulti[leadsNumber];
    float *MIN_PEAK_AMPMulti = ctx->MIN_PEAK_AMPMulti; // static float MIN_PEAK_AMPMulti[leadsNumber];
    int *preBlankCntMulti = ctx->preBlankCntMulti; // static int preBlankCntMulti[leadsNumber];
    int *qpkcntMulti = ctx->qpkcntMulti; // static int qpkcntMulti[leadsNumber];
    int *countMulti = ctx->countMulti; // static int countMulti[leadsNumber];
    int *initBlankMulti = ctx->initBlankMulti; // static int initBlankMulti[leadsNumber];
    float (*rrbufMulti)[8] = ctx->rrbufMulti; //static float rrbufMulti[leadsNumber][8];
    float (*noiseMulti)[8] = ctx->noiseMulti; //static float noiseMulti[leadsNumber][8];
    float *newPeakMulti = ctx->newPeakMulti; // static float newPeakMulti[leadsNumber];
    float (*qrsbufMulti)[8] = ctx->qrsbufMulti; //static float qrsbufMulti[leadsNumber][8];
    float (*rsetBuffMulti)[8] = ctx->rsetBuffMulti; //static float rsetBuffMulti[leadsNumber][8];
    float (*pyMulti)[8] = ctx->pyMulti; //static float pyMulti[leadsNumber][8];
    float (*pxMulti)[8] = ctx->pxMulti; //static float pxMulti[leadsNumber][8];
    float *det_threshMulti = ctx->det_threshMulti; // static float det_threshMulti[leadsNumber];
    int *rsetCountMulti = ctx->rsetCountMulti; // static int rsetCountMulti[leadsNumber];
    float *nmeanMulti = ctx->nmeanMulti; // static float nmeanMulti[leadsNumber];
    float *rrmeanMulti = ctx->rrmeanMulti; // static float rrmeanMulti[leadsNumber];
    float *qmeanMulti = ctx->qmeanMulti; // static float qmeanMulti[leadsNumber];
    int *sblocMulti = ctx->sblocMulti; // static int sblocMulti[leadsNumber];
    int *QrsDelayMulti = ctx->QrsDelayMulti; // static int QrsDelayMulti[leadsNumber];
    float *sbpeakMulti = ctx->sbpeakMulti; // static float sbpeakMulti[leadsNumber];
    float *initMaxMulti = ctx->initMaxMulti; // static float initMaxMulti[leadsNumber];
    int *lastPositionMulti = ctx->lastPositionMulti; // static int lastPositionMulti[leadsNumber];
    int *rMulti = ctx->rMulti; // static int rMulti[leadsNumber];
    int (*qrscountMulti)[maxQueueLength] = ctx->qrscountMulti; //static int qrscountMulti[leadsNumber][maxQueueLength];
    int *adjcentPeakMulti = ctx->adjcentPeakMulti; // static int adjcentPeakMulti[leadsNumber];
    float *currentNewPeakMulti = ctx->currentNewPeakMulti; // static float currentNewPeakMulti[leadsNumber];
    int *lastNewPeakPositionMulti = ctx->lastNewPeakPositionMulti; // static int lastNewPeakPositionMulti[leadsNumber];
    int *currentNewPeakPositionMulti = ctx->currentNewPeakPositionMulti; // static int currentNewPeakPositionMulti[leadsNumber];
    float (*ecgBufferMulti)[maxQueueLength] = ctx->ecgBufferMulti; //static float ecgBufferMulti[leadsNumber][maxQueueLength];
    float *totalBufferMulti = ctx->totalBufferMulti; // static float totalBufferMulti[maxQueueLength];
    float (*PeakPositionMulti)[maxQueueLength] = ctx->PeakPositionMulti; //static float PeakPositionMulti[leadsNumber][maxQueueLength];
    int *MultiVector = ctx->MultiVector; // static int MultiVector[leadsNumber];
    int *multiRRVector = ctx->multiRRVector; // static int multiRRVector[maxQueueLength];

    int        *ctx_RRcount = &(ctx->RRcount);
    int        *ctx_r = &(ctx->r);
    float      *ctx_lostPeakTemp = &(ctx->lostPeakTemp);
    int        *ctx_VectroPointer = &(ctx->VectroPointer);
    int        *ctx_multiRRCount = &(ctx->multiRRCount);
    int        *ctx_sbcount = &(ctx->sbcount);

    int rrValue = 0;
//  11:02  11/02 2023
    if (init) {

        int kkk = 0;
        while (kkk < leadsNumber) {
            indexMulti[kkk] = 0;
            tempPeakMulti[kkk] = 0;
            MovingWindowSumMulti[kkk] = 0;
            MovingWindowSumPtrMulti[kkk] = 1;
            timeSinceMaxMulti[kkk] = 0;
            lastDatumMulti[kkk] = 0;
            maxMulti[kkk] = 0;
            MIN_PEAK_AMPMulti[kkk] = miniPeak;
            preBlankCntMulti[kkk] = 0;
            qpkcntMulti[kkk] = 0;
            countMulti[kkk] = 0;
            initBlankMulti[kkk] = 0;
            newPeakMulti[kkk] = 0;
            det_threshMulti[kkk] = 0;
            rsetCountMulti[kkk] = 0;
            nmeanMulti[kkk] = 0;
            rrmeanMulti[kkk] = 0;
            qmeanMulti[kkk] = 0;
            sblocMulti[kkk] = 0;
            QrsDelayMulti[kkk] = 0;
            sbpeakMulti[kkk] = 0;
            initMaxMulti[kkk] = 0;
            lastPositionMulti[kkk] = 0;
            rMulti[kkk] = 1;
            adjcentPeakMulti[kkk] = 0;
            currentNewPeakMulti[kkk] = 0;
            lastNewPeakPositionMulti[kkk] = 0;
            currentNewPeakPositionMulti[kkk] = 0;
            MultiVector[kkk] = 0; //  % 26/06  2022
            (*ctx_VectroPointer) = 1;
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            int jjj = 0;
            while (jjj < 8) {
                rrbufMulti[kkk][jjj] = MS1000;
                jjj = jjj + 1;
            }// end

            jjj = 0;
            while (jjj < 8) {
                noiseMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            }//

            jjj = 0;
            while (jjj < 8) {
                qrsbufMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            }//

            jjj = 0;
            while (jjj < 8) {
                rsetBuffMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            } // end

            jjj = 0;
            while (jjj < 8) {
                pxMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            }//

            jjj = 0;
            while (jjj < 8) {
                pyMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            } //

            jjj = 0;
            while (jjj < 10) {
                x_dervMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            } //

            jjj = 0;
            while (jjj < 40) {
                dataMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            } //


            jjj = 0;
            while (jjj < maxQueueLength) {
                qrscountMulti[kkk][jjj] = 0;
                ecgBufferMulti[kkk][jjj] = 0;
                PeakPositionMulti[kkk][jjj] = 0;
                jjj = jjj + 1;
            }//

            kkk = kkk + 1;
        } //
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kkk = 0;
        while (kkk < maxQueueLength) {
            totalBufferMulti[kkk] = 0;
            multiRRVector[kkk] = 0;
            kkk = kkk + 1;
        }// end
        (*ctx_multiRRCount) = 8;

        (*ctx_r) = 0;
        (*ctx_RRcount) = 0;
        (*ctx_lostPeakTemp) = 0;
        (*ctx_sbcount) = MS1500;

        return 0;

    }


    indexMulti[multiLeadNumber] = indexMulti[multiLeadNumber] + 1;
    ecgBufferMulti[multiLeadNumber][indexMulti[multiLeadNumber]] = inputdata;
    if (indexMulti[multiLeadNumber] >= (maxQueueLength - 1))                          //   %%  14:04 10/12 2018 ypz
        indexMulti[multiLeadNumber] = 0;
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%[4 30]
    double a[] = {1, -5.48922254603941, 12.6140277755683, -15.5332815784336, 10.8106760640655, -4.03135508957097,
                  0.629160899787621};
    double b[] = {0.0409031291391518, -0.156385624666646, 0.190076246754150, -1.54399425548308e-16, -0.190076246754150,
                  0.156385624666646, -0.0409031291391518};


    pxMulti[multiLeadNumber][0] = inputdata;
    float sum = 0;

    int jj = 0;
    while (jj < 7) {
        sum = sum + (float) b[jj] * pxMulti[multiLeadNumber][jj];
        jj = jj + 1;
    }

    jj = 1;
    while (jj < 7) {
        sum = sum - (float) a[jj] * pyMulti[multiLeadNumber][jj];
        jj = jj + 1;
    }

    int k = 6;
    while (k >= 1) {
        pxMulti[multiLeadNumber][k] = pxMulti[multiLeadNumber][k - 1];
        k = k - 1;
    }

    k = 6;
    while (k >= 2) {
        pyMulti[multiLeadNumber][k] = pyMulti[multiLeadNumber][k - 1];
        k = k - 1;
    }

    pyMulti[multiLeadNumber][1] = sum;
//%    yy(ii)=sum;
    float datum = sum;
    sum = 0;


    int RRInterval = 0;
//% y = 1/8 (2x( nT) + x( nT - T) - x(  nT - 3T) - 2x( nT -4T))
//%  temp1=datum*2+x_derv(4)-x_derv(2)-x_derv(1)*2;
    float temp1 = datum * 2 + x_dervMulti[multiLeadNumber][4] - x_dervMulti[multiLeadNumber][2] -
                  x_dervMulti[multiLeadNumber][1] * 2;
    //for i=1:1:3
    int kk = 0;
    while (kk < 4) {
        x_dervMulti[multiLeadNumber][kk] = x_dervMulti[multiLeadNumber][kk + 1];
        kk = kk + 1;
    }
    x_dervMulti[multiLeadNumber][3] = datum;

    float temp2;
    if (temp1 > 800)
        temp2 = fabsf(temp1) * fabsf(temp1) / 10;
    else if (temp1 > 600)
        temp2 = fabsf(temp1) * fabsf(temp1) / 4;
    else
        temp2 = fabsf(temp1) * fabsf(temp1);

// %%%%%%%%%%%%%%%%%%% Average over an 80 ms Window %%%%%%%%%%%%%%%%%%%
    MovingWindowSumMulti[multiLeadNumber] = MovingWindowSumMulti[multiLeadNumber] + temp2;
    MovingWindowSumMulti[multiLeadNumber] = MovingWindowSumMulti[multiLeadNumber] -
                                            dataMulti[multiLeadNumber][MovingWindowSumPtrMulti[multiLeadNumber]];
    dataMulti[multiLeadNumber][MovingWindowSumPtrMulti[multiLeadNumber]] = temp2;
    MovingWindowSumPtrMulti[multiLeadNumber] = MovingWindowSumPtrMulti[multiLeadNumber] + 1;

    if (MovingWindowSumPtrMulti[multiLeadNumber] > (WINDOW_WIDTH - 1)) {
        MovingWindowSumPtrMulti[multiLeadNumber] = 0;
    }

    float output;
    if ((MovingWindowSumMulti[multiLeadNumber] / WINDOW_WIDTH) > 32000) {
        output = 32000;
    } else {
        output = MovingWindowSumMulti[multiLeadNumber] / WINDOW_WIDTH;
    }
    PeakPositionMulti[multiLeadNumber][indexMulti[multiLeadNumber]] = output;


    float fdatum = output;
    float pk = 0;
    if (timeSinceMaxMulti[multiLeadNumber] > 0) {
        timeSinceMaxMulti[multiLeadNumber] = timeSinceMaxMulti[multiLeadNumber] + 1;
    }

    if ((fdatum > lastDatumMulti[multiLeadNumber]) && (fdatum > maxMulti[multiLeadNumber])) {
        maxMulti[multiLeadNumber] = fdatum;
        if (maxMulti[multiLeadNumber] > 2) {
            timeSinceMaxMulti[multiLeadNumber] = 1;
        }
    } else {
        if (fdatum < maxMulti[multiLeadNumber] / 2) {
            pk = maxMulti[multiLeadNumber];
            maxMulti[multiLeadNumber] = 0;
            timeSinceMaxMulti[multiLeadNumber] = 0;
        } else if (timeSinceMaxMulti[multiLeadNumber] > MS95) {
            pk = maxMulti[multiLeadNumber];
            maxMulti[multiLeadNumber] = 0;
            timeSinceMaxMulti[multiLeadNumber] = 0;
        }
    }
    lastDatumMulti[multiLeadNumber] = fdatum;

    float aPeak = pk; //%% output data
    if (aPeak < MIN_PEAK_AMPMulti[multiLeadNumber]) {
        aPeak = 0;
    }

// %%%%%%%%%%%%%%%%%% peak() end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %% Hold any peak that is detected for 200 ms
    // %% in case a bigger one comes along.  There
    // %% can only be one QRS complex in any 200 ms window.
    newPeakMulti[multiLeadNumber] = 0;
//        if((aPeak!=0)||(index>38))
//            newPeak = 0 ;

    if (aPeak && (preBlankCntMulti[multiLeadNumber] == 0))      //   %%		 If there has been no peak for 200 ms%%
    {
        tempPeakMulti[multiLeadNumber] = aPeak;
        preBlankCntMulti[multiLeadNumber] = PRE_BLANK;     //   %%		// MS200
    } else {
        // if((aPeak==0) &&preBlankCnt)   // %%     // If we have held onto a peak for
        if ((fabsf(aPeak) < EPSILON) &&
            preBlankCntMulti[multiLeadNumber]) {                      //   %%		// 200 ms pass it on  for evaluation.
            //   %%		// 200 ms pass it on  for evaluation.
            preBlankCntMulti[multiLeadNumber] = preBlankCntMulti[multiLeadNumber] - 1;
            if ((preBlankCntMulti[multiLeadNumber]) == 0) {
                newPeakMulti[multiLeadNumber] = tempPeakMulti[multiLeadNumber];
            }//  end
        } else {
            //  if(aPeak)					        // %%	   // If we were holding a peak, but
            if (fabsf(aPeak) > EPSILON) {
                if (aPeak > tempPeakMulti[multiLeadNumber])            // %%	   // start counting to 200 ms again.
                {
                    tempPeakMulti[multiLeadNumber] = aPeak;
                    preBlankCntMulti[multiLeadNumber] = PRE_BLANK;// %%    // MS200//%%    // MS200
                } else {
                    preBlankCntMulti[multiLeadNumber] = preBlankCntMulti[multiLeadNumber] - 1;
                    if ((preBlankCntMulti[multiLeadNumber]) == 0) {
                        newPeakMulti[multiLeadNumber] = tempPeakMulti[multiLeadNumber];
                    }
                    // end
                }// end
            }//
        }// end
    }// end

    if (qpkcntMulti[leadsNumber - 1] < 8) {
        countMulti[multiLeadNumber] = countMulti[multiLeadNumber] + 1;
        if (newPeakMulti[multiLeadNumber] > 0) {
            countMulti[multiLeadNumber] = WINDOW_WIDTH;
        } //
        initBlankMulti[multiLeadNumber] = initBlankMulti[multiLeadNumber] + 1;
        //       %		if(initBlank ==MS720)//%MS1000  %%%%%%%%%%%%%%%%%%%%% ypz
        // if(((initBlank==250)&&((*ctx_RRcount)<1))||initBlank ==MS1000||(((*ctx_RRcount)>=1)&&(initMax>lastPeak/4)&&(initBlank>=(lastPosition-currentPosition1-20))))
        if ((initBlankMulti[multiLeadNumber] >= MS800)) {
            if (initMaxMulti[multiLeadNumber] > MIN_PEAK_AMPMulti[multiLeadNumber]) {
                (*ctx_RRcount) = (*ctx_RRcount) + 1;
            }
            //        initBlank=0;//
            // end
            //      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            float lastInitMax;
            if (((*ctx_RRcount) >= 1))//%&&(newPeak>lastPeak/5))
            {
                //     lastInitMax=lastPeak;
                //     lastPeak=initMax;//%+lastPeak)/((*ctx_RRcount)-2);
                //     if((lastInitMax!=0)&&(lastPeak>lastInitMax*2.5))//%% fix sudden super peak
                //        {
                //         lastPeak=lastPeak/4;
                //          if (lastPeak>lastInitMax*1.5)
                //          {
                //           lastPeak = lastPeak*2;
                //          }//
                //      }//
                int moveCount = 0;
                int movePointer = 0;
                int searchLength = 125;//%78;%46  78 backSearchLength
                int frontSearch = 125;//%72;%39    72
                int currentPosition =
                        lastPositionMulti[multiLeadNumber] - QRSDelayLength - frontSearch;//%-backSearchLength; /
                //    int  currentPosition=index- QRSDelayLength-frontSearch;//%-backSearchLength; /

                //currentPosition = counterAdjust( currentPosition, maxQueueLength );//  %%  14:04 10/12 2018 ypz
                if (currentPosition > maxQueueLength) {
                    currentPosition = currentPosition - maxQueueLength;
                }
                if (currentPosition < 0) {
                    currentPosition = currentPosition + maxQueueLength;
                }
//                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                float currentMax = maxNegative;
                //%     currentPosition=currentPosition-frontSearch;
                int startPosition = currentPosition;
                float goalMax = ecgBufferMulti[multiLeadNumber][currentPosition];
                while (moveCount < searchLength) {
                    if (fabsf(goalMax) > currentMax) {
                        movePointer = moveCount;
                        currentMax = fabsf(goalMax);
                    }

                    currentPosition = currentPosition + 1;
                    if (currentPosition >= maxQueueLength)                    //  %%  14:04 10/12 2018 ypz
                        currentPosition = 0;
                    //
                    goalMax = ecgBufferMulti[multiLeadNumber][currentPosition];
                    moveCount = moveCount + 1;
                }//

                int position = startPosition +
                               movePointer;//= counterAdjust( startPosition+movePointer, maxQueueLength ); // %%  14:04 10/12 2018 ypz

                if (position > maxQueueLength) {
                    position = startPosition + movePointer - maxQueueLength;

                }
                if (position < 0) {
                    position = (startPosition + movePointer) + maxQueueLength;
                }

                qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] = position;//%%startPosition+movePointer;

                MultiVector[(*ctx_VectroPointer)] = 1;
                (*ctx_VectroPointer) = (*ctx_VectroPointer) + 1;
                if ((*ctx_VectroPointer) > leadsNumber) {
                    int ll = 0;
                    while (ll < leadsNumber) {
                        if (MultiVector[ll] == 0)
                            break;
                        ll = ll + 1;
                    }

                    ll = 0;
                    while (ll < leadsNumber) {
                        MultiVector[ll] = 0;
                        ll = ll + 1;
                    }//

                    (*ctx_VectroPointer) = 1;
                }//

//%%%%%%%%%%%%%%%%%%%%%%%%%%  09/07 2022 begin   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                int ll;
                kk = 0;
                while (kk < maxQueueLength) {
                    //  totalPeakPosition[kk]=0;
                    ll = 0;
                    float tempMax = 0;
                    float ecgMax = maxNegative;
                    while (ll < leadsNumber) {
                        if (PeakPositionMulti[ll][kk] > tempMax)
                            tempMax = PeakPositionMulti[ll][kk];

                        if (fabsf(ecgBufferMulti[ll][kk]) > ecgMax)
                            ecgMax = ecgBufferMulti[ll][kk];

                        ll = ll + 1;
                    }
                    // totalPeakPosition(kk)=tempMax;
                    totalBufferMulti[kk] = ecgMax;
                    kk = kk + 1;
                }

                moveCount = 0;
                movePointer = 0;
                searchLength = 125;//%78;%46  78 backSearchLength
                frontSearch = 125;//%72;%39    72
                currentPosition =
                        lastPositionMulti[multiLeadNumber] - QRSDelayLength - frontSearch;//%-backSearchLength;
                //% currentPosition=indexMulti(multiLeadNumber)- QRSDelayLength-frontSearch;//%-backSearchLength
                currentPosition = counterAdjust(currentPosition,
                                                maxQueueLength); // %%  14:04 10/12 2018 ypz
                currentMax = maxNegative;
                //%     currentPosition=currentPosition-frontSearch;
                startPosition = currentPosition;
                // %% goalMax=totalBufferMulti(multiLeadNumber,currentPosition);
                goalMax = totalBufferMulti[currentPosition];
                while (moveCount < searchLength) {
                    if (fabsf(goalMax) > currentMax) {
                        movePointer = moveCount;
                        currentMax = fabsf(goalMax);
                    }
                    currentPosition = currentPosition + 1;
                    if (currentPosition == maxQueueLength + 1)                    //  %%  14:04 10/12 2018 ypz
                        currentPosition = 1;
                    //%% goalMax=totalBufferMulti(multiLeadNumber,currentPosition);
                    goalMax = totalBufferMulti[currentPosition];
                    moveCount = moveCount + 1;
                }

                position = counterAdjust(startPosition + movePointer, maxQueueLength);

//%%%%%%%%%%%%%%%%%%%%%%%%  09/07 2022  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                adjcentPeakMulti[multiLeadNumber] = indexMulti[multiLeadNumber];
//         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 09:30 05/03 2018

                if (RRInterval > maxQueueLength) //  %% circular queue                        %%  14:04 10/12 2018 ypz
                {
                    RRInterval = 0;
                }
                // %    RRInterval =index-RRInterval;
                if (RRInterval < 0) {
                    RRInterval = RRInterval + maxQueueLength;
                }
                rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + 1;
                if (rMulti[multiLeadNumber] ==
                    maxQueueLength + 1) {  // %% circular queue                        %%  14:04 10/12 2018 ypz
                    rMulti[multiLeadNumber] = 1;
                }
            }

//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            //  15:11 17/08 2017 ypz
            initBlankMulti[multiLeadNumber] = 0;
            if (initMaxMulti[multiLeadNumber] > 20000) {
                qrsbufMulti[multiLeadNumber][qpkcntMulti[multiLeadNumber]] = initMaxMulti[multiLeadNumber] / 3;
            } else {
                qrsbufMulti[multiLeadNumber][qpkcntMulti[multiLeadNumber]] = initMaxMulti[multiLeadNumber];
            }
            initMaxMulti[multiLeadNumber] = 0;
            qpkcntMulti[multiLeadNumber] = qpkcntMulti[multiLeadNumber] + 1;//%qpkcnt=qpkcnt+1 ;

//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (qpkcntMulti[leadsNumber - 1] == 8) {
                nmeanMulti[multiLeadNumber - 1] = 0;
                jj = 1;
                while (jj < leadsNumber) {
                    rrmeanMulti[jj] = MS1000;
                    //sbcountMulti[jj]= MS1500+MS150 ;
                    jj = jj + 1;
                }

                //det_threshMulti= det_threshMultiCompute(qrsbufMulti,nmeanMulti,TH ,7,leadsNumber);
                int ii = 0;
                while (ii < leadsNumber) {
                    jj = 0;
                    sum = 0;
                    while (jj < 8) {
                        sum = sum + qrsbufMulti[ii][jj];
                        jj = jj + 1;
                    }
                    qmeanMulti[ii] = sum / 8;
                    float temp = qmeanMulti[ii] - nmeanMulti[ii];
                    temp = temp * (float) TH;
                    float dmed = temp;
                    float thrsh = nmeanMulti[ii] + dmed;
                    det_threshMulti[ii] =
                            WeightFactor * qmeanMulti[ii] + (1 - ForgetFactor) * thrsh; //%adaptive threshold
                    ii = ii + 1;
                }
            }

        }


        if (newPeakMulti[multiLeadNumber] > initMaxMulti[multiLeadNumber]) {
            lastPositionMulti[multiLeadNumber] = indexMulti[multiLeadNumber];
            initMaxMulti[multiLeadNumber] = newPeakMulti[multiLeadNumber];
        }

    } else {
//  > 8 beats
        countMulti[multiLeadNumber] = countMulti[multiLeadNumber] + 1;

        if (newPeakMulti[multiLeadNumber] > 0) {
            int blankNumber = indexMulti[multiLeadNumber] - adjcentPeakMulti[multiLeadNumber];//%%adjcentPeak;
            if (blankNumber < 0)
                blankNumber = blankNumber + ModNumber;
            //
            if ((newPeakMulti[multiLeadNumber] > det_threshMulti[multiLeadNumber]) && (blankNumber > MS200))
                //	{
            {
                adjcentPeakMulti[multiLeadNumber] = indexMulti[multiLeadNumber];
                currentNewPeakMulti[multiLeadNumber] = newPeakMulti[multiLeadNumber];
                lastNewPeakPositionMulti[multiLeadNumber] = currentNewPeakPositionMulti[multiLeadNumber];
                currentNewPeakPositionMulti[multiLeadNumber] = indexMulti[multiLeadNumber];
                //%%%memmove(&qrsbuf[1], qrsbuf, MEMMOVELEN) ;
                float qrsTemp = 0;
                int ii = 0;
                while (ii < 7) {
                    qrsTemp = qrsTemp + qrsbufMulti[multiLeadNumber][ii];
                    qrsbufMulti[multiLeadNumber][ii] = qrsbufMulti[multiLeadNumber][ii + 1];
                    ii++;
                }//

                qrsbufMulti[multiLeadNumber][7] = newPeakMulti[multiLeadNumber];

                sum = 0;
                ii = 0;
                while (ii < 7) {
                    sum = sum + qrsbufMulti[multiLeadNumber][ii];
                    ii = ii + 1;
                }//
                sum = sum / 8;
                qmeanMulti[multiLeadNumber] = sum;

                float thrsh = 0;
                float dmed = 0;
                float temp = 0;
                dmed = qmeanMulti[multiLeadNumber] - nmeanMulti[multiLeadNumber];
                temp = dmed;
                temp = temp * (float) TH;
                dmed = temp;
                thrsh = nmeanMulti[multiLeadNumber] + dmed;
                //%det_thresh=thrsh ;%WeightFactor
                det_threshMulti[multiLeadNumber] =
                        WeightFactor * qmeanMulti[multiLeadNumber] + (1 - ForgetFactor) * thrsh;

                //%%%memmove(&rrbuf[1], rrbuf, MEMMOVELEN) ;
                ii = 0;
                while (ii < 7) {
                    rrbufMulti[multiLeadNumber][ii] = rrbufMulti[multiLeadNumber][ii + 1];
                    ii = ii + 1;
                }//
                rrbufMulti[multiLeadNumber][7] = (float) countMulti[multiLeadNumber] - WINDOW_WIDTH;

                //%%%rrmean = mean(rrbuf,8) ;
                ii = 0;
                while (ii < 7) {
                    sum = sum + rrbufMulti[multiLeadNumber][ii];
                    ii = ii + 1;
                }//
                sum = sum / 8;
                rrmeanMulti[multiLeadNumber] = sum;

                (*ctx_sbcount) = (int) (rrmeanMulti[multiLeadNumber] + (rrmeanMulti[multiLeadNumber] / 2) + WINDOW_WIDTH);
                countMulti[multiLeadNumber] = WINDOW_WIDTH;
                sbpeakMulti[multiLeadNumber] = 0;
                QrsDelayMulti[multiLeadNumber] = WINDOW_WIDTH + FILTER_DELAY;

//%%%%%%%%%%%%%%%%%%%%%%%% QRS positin adjust %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                int moveCount = 0;
                int movePointer = 0;
                int searchLength = 120;//%+140;%120;%100;%120;%36
                int frontSearch = 115;//%80;%72;%30

                int currentPosition = indexMulti[multiLeadNumber] - QRSDelayLength - frontSearch;
                if (currentPosition < 1)
                    currentPosition = currentPosition + maxQueueLength;
                //
                float currentMax = maxNegative;//%-32768;
                float startPosition = (float) currentPosition;
                float goalMax = (ecgBufferMulti[multiLeadNumber][currentPosition]);//%%abs(y(currentPosition));  BeatBuffer global ECGBuffer;
                float maxQRS = (ecgBufferMulti[multiLeadNumber][currentPosition]);
                float minQRS = (ecgBufferMulti[multiLeadNumber][currentPosition]);
                while (moveCount < searchLength) {
                    if (fabsf(goalMax) > currentMax) {
                        movePointer = moveCount;
                        currentMax = fabsf(goalMax);
                    }
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if (ecgBufferMulti[multiLeadNumber][currentPosition] > maxQRS)
                        maxQRS = ecgBufferMulti[multiLeadNumber][currentPosition];
                    else if (ecgBufferMulti[multiLeadNumber][currentPosition] < minQRS)
                        minQRS = ecgBufferMulti[multiLeadNumber][currentPosition];

                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    currentPosition = currentPosition + 1;
                    if (currentPosition == ModNumber + 1)
                        currentPosition = 1;

                    goalMax = ecgBufferMulti[multiLeadNumber][currentPosition];
                    moveCount = moveCount + 1;
                }

                RRInterval = (int) (startPosition + (float) movePointer);
                if (RRInterval > ModNumber)
                    RRInterval = RRInterval - ModNumber;
                else if (RRInterval == ModNumber + 1)
                    RRInterval = 1;

                // ooo=multiLeadNumber;
                //%               qrscountMulti(multiLeadNumber,rMulti(multiLeadNumber))=RRInterval;%startPosition+movePointer;
                qrscountMulti[multiLeadNumber][(*ctx_multiRRCount)] = RRInterval;//%startPosition+movePointer;
                float Crit = 0.02f * (maxQRS - minQRS);


                if (rMulti[multiLeadNumber] > 1) {
                    int outNumber = rMulti[multiLeadNumber] - 1;
                    if (outNumber < 1)
                        outNumber = outNumber + maxQueueLength;

                    // %%rrValue=qrscount((*ctx_r)-1);  %%%%%%%%%%%%%%10:20 23/04 2018
                    //%            rrValue=qrscountMulti(multiLeadNumber,outNumber);
                }

                initBlankMulti[multiLeadNumber] = 0;
                initMaxMulti[multiLeadNumber] = 0;
                rsetCountMulti[multiLeadNumber] = 1;
                //%               if(noiseFlag==0)
                rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + 1;

                if (rMulti[multiLeadNumber] == ModNumber + 1)
                    rMulti[multiLeadNumber] = 1;

                // %%%%%%%%%%%%%%%%%%%%%%%%%%  09/07 2022 begin   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (multiLeadNumber >= (leadsNumber - 1)) {
                    kk = 0;
                    while (kk < maxQueueLength) {
                        //float totalPeakPosition[kk]=0;
                        //int ll=1;
                        float peakMax = maxNegative;
                        float peakMin = - maxNegative;
                        float ecgMax = maxNegative;
                        // % ecgMin=-maxNegative;

                        int ll = 0;
                        float peakSum = 0;
                        float ecgSum = 0;
                        while (ll < leadsNumber) {
                            if (PeakPositionMulti[ll][kk] > peakMax)
                                peakMax = PeakPositionMulti[ll][kk];
                            else if (PeakPositionMulti[ll][kk] < peakMin)
                                peakMin = PeakPositionMulti[ll][kk];
                            else
                                peakSum = peakSum + PeakPositionMulti[ll][kk];

                            ecgSum = ecgSum + fabsf(ecgBufferMulti[ll][kk]);

                            ll = ll + 1;
                        }

                        totalBufferMulti[kk] = ecgSum / (leadsNumber - 4);
                        kk = kk + 1;
                    }

                    moveCount = 0;
                    movePointer = 0;
                    searchLength = 125;//%78;%46  78 backSearchLength
                    frontSearch = 125;//%72;%39    72
                    currentPosition = indexMulti[multiLeadNumber] - QRSDelayLength - frontSearch;//%-backSearchLength;
                    // % currentPosition=indexMulti(multiLeadNumber)- QRSDelayLength-frontSearch;%-backSearchLength
                    currentPosition = counterAdjust(currentPosition,
                                                    maxQueueLength);  //%%  14:04 10/12 2018 ypz
                    currentMax = maxNegative;
                    // %     currentPosition=currentPosition-frontSearch;
                    startPosition = (float) currentPosition;
                    // %% goalMax=totalBufferMulti(multiLeadNumber,currentPosition);
                    goalMax = totalBufferMulti[currentPosition];
                    while (moveCount < searchLength) {
                        if (fabsf(goalMax) > currentMax) {
                            movePointer = moveCount;
                            currentMax = fabsf(goalMax);
                        }
                        currentPosition = currentPosition + 1;
                        if (currentPosition == maxQueueLength + 1)                   //   %%  14:04 10/12 2018 ypz
                            currentPosition = 0;

                        // %% goalMax=totalBufferMulti(multiLeadNumber,currentPosition);
                        goalMax = totalBufferMulti[currentPosition];
                        moveCount = moveCount + 1;
                    }//

                    float positionpp = (float) counterAdjust((int) (startPosition + (float) movePointer),
                                                             maxQueueLength);

                    multiRRVector[(*ctx_multiRRCount)] = (int) positionpp;
                    temp = (float) indexMulti[multiLeadNumber] - positionpp;
                    if (temp < 0)
                        temp = temp + maxQueueLength;
                    //
                    //%         differenceFrequencyTime((*ctx_multiRRCount))=temp;
                    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    //%     if((*ctx_multiRRCount)>1)
                    kk = 0;
                    int kkCount = 0;
                    int IntervalDiff = 20;//%%  20
                    int pp = ((*ctx_multiRRCount) - 1);
                    if ((pp == maxQueueLength) || (pp < 0))                    //  %%  14:04 10/12 2018 ypz
                        pp = 0;
                    //
                    int CurrentInterval = multiRRVector[pp];
                    temp1 = (float) (CurrentInterval - IntervalDiff);
                    if (temp1 < 0)
                        temp1 = temp1 + maxQueueLength;
                    //
                    temp2 = (float) (CurrentInterval + IntervalDiff);
                    if (temp1 > maxQueueLength)
                        temp2 = temp2 - maxQueueLength;
                    //
                    while ((pp >= 8) && (kk < leadsNumber)) {

                        if ((((float) qrscountMulti[kk][pp] > temp1) && ((float) qrscountMulti[kk][pp] < temp2)) ||
                            (((float) qrscountMulti[kk][pp + 1] >= temp1) &&
                             ((float) qrscountMulti[kk][pp + 1] <= temp2)))
                            kkCount = kkCount + 1;
                        //
                        kk = kk + 1;
                    }//

                    if ((kkCount > 2) || (CurrentInterval < (IntervalDiff + 2)) ||
                        (CurrentInterval > (maxQueueLength - IntervalDiff)))
                        rrValue = CurrentInterval;
                    else {
                        if (qrscountMulti[11][pp] > 0) {
                            rrValue = qrscountMulti[11][pp];
                        } else if (qrscountMulti[2][pp] > 0) {
                            rrValue = qrscountMulti[2][pp];
                        }//
                    }//

                    (*ctx_multiRRCount) = (*ctx_multiRRCount) + 1;
                    if ((*ctx_multiRRCount) == maxQueueLength + 1) {                   //   %%  14:04 10/12 2018 ypz
                        (*ctx_multiRRCount) = 1;
                    }//

                }// %%   end if(multiLeadNumber==leadsNumber)
            }
//  $$$$$$$$$$$$$$$$$$$$$$$$$ for noise update $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            else  //  second else  //noise
            {
                //%%% memmove(&noise[1],noise,MEMMOVELEN) ;
                float noiseTemp = 0;
                for (int jj = 0; jj < 7; jj++) {
                    noiseTemp = noiseTemp + noiseMulti[multiLeadNumber][jj];
                    noiseMulti[multiLeadNumber][jj] = noiseMulti[multiLeadNumber][jj + 1];//% noise(i+1);
                }//  end
                if (newPeakMulti[multiLeadNumber] > 2000) {
                    noiseTemp = noiseTemp / 7;
                    noiseTemp = (noiseTemp * 17 + newPeakMulti[multiLeadNumber]) / 18;
                    noiseMulti[multiLeadNumber][7] = noiseTemp;
                } else {
                    noiseMulti[multiLeadNumber][7] = newPeakMulti[multiLeadNumber];
                }//
                //%%%	nmean = mean(noise,8) ;
                sum = 0;
                for (int jj = 0; jj < 7; jj++) {
                    sum = sum + noiseMulti[multiLeadNumber][jj];
                }//
                sum = sum / 8;
                nmeanMulti[multiLeadNumber] = sum;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                //%%%				det_thresh = thresh(qmean,nmean) ;
                float thrsh = 0;
                float dmed = 0;
                float temp = 0;
                dmed = qmeanMulti[multiLeadNumber] - nmeanMulti[multiLeadNumber];
                temp = dmed;
                temp = temp * (float) TH;
                dmed = temp;
                thrsh = nmeanMulti[multiLeadNumber] + dmed;
                det_threshMulti[multiLeadNumber] = thrsh * 1.2f; //%%%%%%%%%%%%%%%%%%%%%%% 09:31 2018
                //%%%					// Don't include early peaks (which might be T-waves)
                //%%%					// in the search back process.  A T-wave can mask
                //%%%					// a small following QRS.
                if ((newPeakMulti[multiLeadNumber] > sbpeakMulti[multiLeadNumber]) &&
                    ((countMulti[multiLeadNumber] - WINDOW_WIDTH) >= MS360)) {
                    sbpeakMulti[multiLeadNumber] = newPeakMulti[multiLeadNumber];
                    sblocMulti[multiLeadNumber] = countMulti[multiLeadNumber] - WINDOW_WIDTH;
                    //lostPeakTempMulti[multiLeadNumber]=indexMulti[multiLeadNumber];
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                }//   if((newPeakMulti(multiLeadNumber) > sbpeakMulti(multiLeadNumber)) && ((countMulti(multiLeadNumber)-WINDOW_WIDTH) >= MS360))

            }// end else noise

            //		}// end  if((newPeakMulti(multiLeadNumber)
        }// end if(newPeakMulti[multiLeadNumber]>0)


        if ((rMulti[multiLeadNumber] > 2) && (countMulti[multiLeadNumber] > (*ctx_sbcount)) &&
            (det_threshMulti[multiLeadNumber] > MIN_PEAK_AMP) && (sbpeakMulti[multiLeadNumber] >
                                                                  (det_threshMulti[multiLeadNumber] /
                                                                   2)))  //%%  0.4  0.23   13:31 13/06 2018 ypz
        {
//%if((count > (*ctx_sbcount)*0.4) && (sbpeak > (det_thresh*0.23)))%%  0.4  0.23   13:31 13/06 2018 ypz
            // temp=countMulti[multiLeadNumber];
//%if((count > (*ctx_sbcount)*0.88) && (sbpeak > (det_thresh*0.23)))
            for (int jj = 0; jj < 7; jj++) {
                qrsbufMulti[multiLeadNumber][jj] = qrsbufMulti[multiLeadNumber][jj + 1];
            }//
            qrsbufMulti[multiLeadNumber][7] = sbpeakMulti[multiLeadNumber];//%newPeak ;
            //%%%			qmean = mean(qrsbuf,8) ;
            sum = 0;
            for (int ii = 0; ii < 8; ii++) {
                sum = sum + qrsbufMulti[multiLeadNumber][ii];
            }//
            sum = sum / 8;
            qmeanMulti[multiLeadNumber] = sum;
            //%%%			det_thresh = thresh(qmean,nmean) ;
            float thrsh = 0;
            float dmed = 0;
            float temp = 0;
            dmed = qmeanMulti[multiLeadNumber] - nmeanMulti[multiLeadNumber];
            temp = dmed;
            temp = temp * (float) TH;
            dmed = temp;
            thrsh = nmeanMulti[multiLeadNumber] + dmed;
            det_threshMulti[multiLeadNumber] = thrsh;
            //%%%				memmove(&rrbuf[1],rrbuf,MEMMOVELEN) ;
            for (int jj = 0; jj < 7; jj++) {
                rrbufMulti[multiLeadNumber][jj] = rrbufMulti[multiLeadNumber][jj + 1];
            }//
            rrbufMulti[multiLeadNumber][7] = (float) sblocMulti[multiLeadNumber];
            //%%%	     rrbuf(1)= sbloc ;
            //%%%		 rrmean = mean(rrbuf,8) ;
            sum = 0;
            for (int jj = 0; jj < 8; jj++) {
                sum = sum + rrbufMulti[multiLeadNumber][jj];
            }// end
            sum = sum / 8;
            rrmeanMulti[multiLeadNumber] = sum;
            (*ctx_sbcount) = (int) (rrmeanMulti[multiLeadNumber] + rrmeanMulti[multiLeadNumber] / 2 + (float) WINDOW_WIDTH);
            QrsDelayMulti[multiLeadNumber] = countMulti[multiLeadNumber] - sblocMulti[multiLeadNumber];
            countMulti[multiLeadNumber] = countMulti[multiLeadNumber] - sblocMulti[multiLeadNumber];
            // if(QrsDelayMulti[multiLeadNumber]>3000)
            //    kkkk=0;
            //
            if (QrsDelayMulti[multiLeadNumber] > maxQueueLength)
                // %             rrValue=maxQueueLength-1;
                return 0;
            //
            QrsDelayMulti[multiLeadNumber] = QrsDelayMulti[multiLeadNumber] + FILTER_DELAY;

            //  %   SearchBackRRcount=SearchBackRRcount+1;
            //%         moveCount=0;
            float moveCount = 0;
            float movePointer = 0;
            int searchLength = 110;//%36
            int frontSearch = 110;//%30
            // % currentPosition=indexMuiti(muliLeadNumber)- QRSDelayLength-frontSearch+9;
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            int currentPosition = indexMulti[multiLeadNumber] - QrsDelayMulti[multiLeadNumber] - frontSearch;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (currentPosition < 1)
                currentPosition = currentPosition + maxQueueLength;
            //
            float currentMax = -32768;
            //%    currentPosition=currentPosition-frontSearch;
            int startPosition = currentPosition;
            //if(currentPosition>4000)
            //OOO=0;
            //
            float goalMax = fabsf(
                    ecgBufferMulti[multiLeadNumber][currentPosition]);//%%abs(y(currentPosition));  BeatBuffer global ECGBuffer;

            while (moveCount < (float) searchLength) {
                if (goalMax > currentMax) {
                    movePointer = moveCount;
                    currentMax = goalMax;
                }//
                currentPosition = currentPosition + 1;
                if (currentPosition == maxQueueLength + 1)
                    currentPosition = 1;
                //
                goalMax = fabsf(ecgBufferMulti[multiLeadNumber][currentPosition]);
                moveCount = moveCount + 1;
            }//
            RRInterval = (int) ((float) startPosition + movePointer);
            if (RRInterval > ModNumber) {
                RRInterval = RRInterval - ModNumber;
            }//
            //%  qrscountMulti(multiLeadNumber,rMulti(multiLeadNumber))=RRInterval;

            MultiVector[(*ctx_VectroPointer)] = 1;
            (*ctx_VectroPointer) = (*ctx_VectroPointer) + 1;
            if ((*ctx_VectroPointer) >= leadsNumber) {
                int ll = 0;
                while (ll < leadsNumber) {
                    MultiVector[(*ctx_VectroPointer)] = 0;
                    ll = ll + 1;
                }//
                (*ctx_VectroPointer) = 1;
            }//


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            /*
                 xx=zeros(1,2*MS200+MS280);
                 jjStart=RRInterval-2*MS200;
                 if(jjStart<1)
                     jjStart=jjStart+maxQueueLength;
                 end
                 jjEnd=RRInterval+MS280;
                 if(jjEnd>maxQueueLength)
                    jjEnd=jjEnd-maxQueueLength;
                 end

                xx= ecgCopy(ECGBuffer,jjStart,jjEnd);
                noiseLevel= HFNoiseCheck(xx);
                if(noiseLevel>75)
                  qrscountMulti(multiLeadNumber,rMulti(multiLeadNumber))=0;
                end
					*/
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            float qrsDiff = (float) (qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] -
                                     qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber] - 1]);
            if (qrsDiff < 0)
                qrsDiff = qrsDiff + ModNumber;
            //

            //%     if (qrscount((*ctx_r))-qrscount((*ctx_r)-1))<125
            if (qrsDiff < 125) {
                int qrsPre = rMulti[multiLeadNumber] - 1;
                if (qrsPre < 0)
                    qrsPre = qrsPre + ModNumber;
                //  end
                // %       if qrscount((*ctx_r)) < qrscount((*ctx_r)-1)
                //%      if qrscountMulti(multiLeadNumber,rMulti(multiLeadNumber)) < qrscount(qrsPre)
                if (qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] < qrscountMulti[multiLeadNumber][qrsPre]) {
                    float qrsValue = (float) qrscountMulti[multiLeadNumber][qrsPre];
                    rMulti[multiLeadNumber] = rMulti[multiLeadNumber] - 1;
                    if (rMulti[multiLeadNumber] < 0) {
                        rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + ModNumber;
                    }//
                    qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] = (int) qrsValue;
                } else if ((qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber] - 1] != 0) &&
                           (qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] != 0) &&
                           (fabsf(ecgBufferMulti[multiLeadNumber][qrscountMulti[multiLeadNumber][rMulti[
                                   multiLeadNumber -
                                   1]]]) <
                            fabsf(ecgBufferMulti[multiLeadNumber][qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]]]))) {
                    float qrsValue = (float) qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]];
                    rMulti[multiLeadNumber] = rMulti[multiLeadNumber] - 1;
                    if (rMulti[multiLeadNumber] < 0) {
                        rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + ModNumber;
                    }//
                    qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] = (int) qrsValue;
                } else {
                    qrscountMulti[multiLeadNumber][rMulti[multiLeadNumber]] = 0;
                    rMulti[multiLeadNumber] = rMulti[multiLeadNumber] - 1;
                    if (rMulti[multiLeadNumber] < 0) {
                        rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + ModNumber;
                    }//
                }//
            }//   if(qrsDiff<125)

            if (rMulti[multiLeadNumber] > 0) {
                int outNumber = rMulti[multiLeadNumber] - 1;
                if (outNumber < 0) {
                    // TODO: It's suspected to be a bug, waiting for repair.
                    outNumber = outNumber + maxQueueLength;
                }//
                //%              rrValue=qrscountMulti(multiLeadNumber,outNumber);
            }// end
            if (RRInterval == maxQueueLength + 1)  // %% circular queue
            {
                RRInterval = 1;
            }//
            RRInterval = indexMulti[multiLeadNumber] - RRInterval;
            if (RRInterval < 1) {
                RRInterval = RRInterval + maxQueueLength;
            }//

            //%    monitoringPointerMulti(multiLeadNumber)=0;
            // % monintorArray=zeros(1,6);

            // %lastPeakIndex=adjcentPeak;
            //adjcentPeak=(*ctx_lostPeakTemp);
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rMulti[multiLeadNumber] = rMulti[multiLeadNumber] + 1;
            // qrsFlag=1;
            sbpeakMulti[multiLeadNumber] = 0;
            initBlankMulti[multiLeadNumber] = 0;
            initMaxMulti[multiLeadNumber] = 0;
            rsetCountMulti[multiLeadNumber] = 1;

        }// if((rMulti(multiLeadNumber)>2)&&(countMulti(multiLeadNumber)>(*ctx_sbcount))&&(det_threshMulti(multiLeadNumber)>MIN_PEAK_AMP)&& (sbpeakMulti(multiLeadNumber)>(

    }// Big else if(qpkcntMulti[leadsNumber]< 8)

    if (qpkcntMulti[multiLeadNumber] == 8) {
        initBlankMulti[multiLeadNumber] = initBlankMulti[multiLeadNumber] + 1;
        if (initBlankMulti[multiLeadNumber] == MS1000) {
            initBlankMulti[multiLeadNumber] = 0;
            rsetBuffMulti[multiLeadNumber][rsetCountMulti[multiLeadNumber]] = initMaxMulti[multiLeadNumber];
            initMaxMulti[multiLeadNumber] = 0;
            rsetCountMulti[multiLeadNumber] = rsetCountMulti[multiLeadNumber] + 1;

            if (rsetCountMulti[multiLeadNumber] == 8) {
                for (int i = 0; i < 8; i++)//for i = 1:1: 8
                {
                    qrsbufMulti[multiLeadNumber][i] = rsetBuffMulti[multiLeadNumber][i];
                    noiseMulti[multiLeadNumber][i] = 0;
                }// end
                //%%%				qmean = mean( rsetBuff, 8 ) ;
                sum = 0;
                for (int j = 0; j < 8; j++)//for i=1:1:8
                {
                    sum = sum + rsetBuffMulti[multiLeadNumber][j];
                }//
                sum = sum / 8;
                qmeanMulti[multiLeadNumber] = sum;
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nmeanMulti[multiLeadNumber] = 0;
                rrmeanMulti[multiLeadNumber] = MS1000;
                (*ctx_sbcount) = MS1500 + MS150;
                //%%%	det_thresh = thresh(qmean,nmean) ;
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                float thrsh = 0;
                float dmed = 0;
                float temp = 0;
                dmed = qmeanMulti[multiLeadNumber] - nmeanMulti[multiLeadNumber];
                temp = dmed;
                temp = temp * (float) TH;
                dmed = temp;
                thrsh = nmeanMulti[multiLeadNumber] + dmed;
                det_threshMulti[multiLeadNumber] = thrsh;
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                initBlankMulti[multiLeadNumber] = 0;
                initMaxMulti[multiLeadNumber] = 0;
                rsetCountMulti[multiLeadNumber] = 1;

            }// end  if(rsetCountMulti[multiLeadNumber]== 8)

        }//	 if(initBlankMulti[multiLeadNumber]== MS1000)

        if (newPeakMulti[multiLeadNumber] > initMaxMulti[multiLeadNumber]) {
            initMaxMulti[multiLeadNumber] = newPeakMulti[multiLeadNumber];
        }


    }//  if( qpkcntMulti[multiLeadNumber]== 8 )

    return rrValue;

} // int qrsDetect2023(float inputdata,int multiLeadNumber,float*ECGBuffer,int init)
		
		
		
		
		
		
		
		
