#include "utility.h"
#include <math.h>

#define  maxQueueLength 3000
#define  BEAT_MS10 5
#define  BEAT_MS50 25
#define  AVELENGTH BEAT_MS50
#define  BEAT_MS70  35
#define  BEAT_MS80  40
#define  BEAT_MS280 140
#define  dataLength 500
#define  BEAT_MS50  25
#define  BEAT_MS110 55
#define  FIDMARK    200


//void utility::LowPassFilter(const float*X,float*Xpc)
void LowPassFilter(const float *X, float *Xpc) {
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


// int utility::IsoCheckDiffRight(float*data,int j,int isoLength,float limit )
int IsoCheckDiffRight(float *data, int j, int isoLength, float limit) {
    // %   IsoCheck determines whether the amplitudes of a run
    // %   of data fall within a sufficiently sall amplitude that
    // %	the run can be considered isoelectric.

    int i = j;
    int k = 0;
    // int returnValue=0;
    //       float  max=data[i];
    //       float  min=data[i];
    while (k < isoLength) {
        //   %  if(abs(data(i-k)-data(i-k+1))>=limit)
        if (fabsf(data[i + k] - data[i + k - 1]) >= limit / 2)
            // returnValue=1;
            return 1;
        // end

//              if(data[i+k]>max)
//                  max=data[i+k];
//              else
//                 if(data[i+k]<min)
//                    min=data[i+k];
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


// float utility:: ECGMin(float*array,int length )
float ECGMin(const float *array, int length) {
    float returnMin;
    if (length >= 1) {
        returnMin = array[0];
        int kk = 0;
        while (kk < length) {
            if (array[kk] < returnMin)
                returnMin = array[kk];
            // end
            kk = kk + 1;
        }//end
    } else {
        returnMin = 0;
    }// end

    return returnMin;

}


//float utility:: ecgSum(float*array,int len)
float ecgSum(float *array, int len) {


    int index = 1;
    float artifactsIndex = 0;
    while (index <= len) {
        artifactsIndex = artifactsIndex + fabsf(array[index]);
        index = index + 1;
    }//end
    artifactsIndex = artifactsIndex / (float) len;

    return artifactsIndex;


}


//float utility::ECGMax(float*array,int length )
float ECGMax(const float *array, int length) {
    float returnMax;
    if (length >= 1) {
        returnMax = array[0];
        int kk = 1;
        while (kk < length) {
            if (array[kk] > returnMax)
                returnMax = array[kk];
            //end
            kk = kk + 1;
        }//end
    } else
        returnMax = 0;
    // end

    return returnMax;
}


//   /////////////////////////////////////////////////////////////////////
// void utility::ecgCopy(float*array,float*result,int start,int finish)
void ecgCopy(const float *array, float *result, int start, int finish) {
    //maxQueueLength=4000;
    int len = finish - start;
    if (len < 0)
        len = len + maxQueueLength;
    // end
    int ii = 0;
    while (ii < len) {
        result[ii] = 0;
        ii = ii + 1;
    }
    //ecgArray=zeros(1,len);

    int jj = 0;
    int index = start;
    while (jj < len) {
        result[jj] = array[index];
        jj = jj + 1;
        index = index + 1;

        if (index > (maxQueueLength - 1))
            index = 0;
        // end
    }

    // end

    // end

}

//   //////////////////////////////////////////////////////
float ecgMean(const float *array, int len) {

    float ecgSum = 0;
    int index = 0;
    while (index <= len) {
        ecgSum = ecgSum + array[index];
        index = index + 1;

    }//end
    return ecgSum / (float) len;

    // end



}

//   /////////////////////////////////////////////////////////////
void ecgDiff(const float *array, float *result, int len) {

    // len=length(array);
    //  diffArray=zeros(1,len);
    int ii = 0;
    while (ii < len) {
        result[ii] = 0;
        ii++;
    }

    int index = 0;
    while (index < (len - 1)) {
        result[index] = array[index + 1] - array[index];
        index = index + 1;
    }
    //   end
    //   end


}

//   /////////////////////////////////////////////////
float noiseLevel(const float *array, int len) {

    int ifinal;
    int iew = len;//length(X);
    float pNoise = 0;
    int jj = 0;
    int ii = 0;
    int kk = 0;
    float ymax = 0;
    float ymin = 0;
    while (jj < iew) {
        ii = ii + 1;
        if ((jj + BEAT_MS10) < iew)
            ifinal = (jj + BEAT_MS10);
        else
            ifinal = iew;
        // end
        ymax = array[jj];
        ymin = array[jj];
        kk = jj;//%ifinal;
        while ((kk < iew) && (kk <= ifinal)) {
            if (array[kk] < ymin)
                ymin = array[kk];
            // end
            if (array[kk] > ymax)
                ymax = array[kk];
            // end
            kk = kk + 1;
        }// end
        pNoise = (ymax - ymin) + pNoise;
        jj = ifinal;
    }// end   while((kk<iew)&&(kk<=ifinal))
    if (ii > 0)
        pNoise = fabsf(pNoise / (float) ii);
    // end
    return pNoise;
}

//   /////////////////////////////////////////////////////


float HFnoiseCheck(float *array) {

    float maxNoiseAve = 0;
    float sum = 0;
    int avePtr = 0;
    float qrsMax = 0;
    float qrsMin = 0;
    float returnValue;
    float aveBuff[AVELENGTH];
    //aveBuff=zeros(1,AVELENGTH);
    int kk = 0;
    while (kk < AVELENGTH) {
        aveBuff[kk] = 0;
        kk = kk + 1;
    }

    float ppArray[500] = {0};
    int ii = 0;
    while (ii < 399) {
        ppArray[ii] = array[ii];
        ii++;
    }




    // %% Determine the QRS amplitude.
    int i = FIDMARK - BEAT_MS70;
    int limit = (FIDMARK + BEAT_MS80);
    while (i < limit) {
        if (array[i] > qrsMax)
            qrsMax = array[i];
        else if (array[i] < qrsMin)
            qrsMin = array[i];
        //end
        // end
        i = i + 1;
    }//end

    i = FIDMARK - BEAT_MS280;
    while (i < FIDMARK + BEAT_MS280) {
        sum = sum - aveBuff[avePtr];
        aveBuff[avePtr] = fabsf(array[i] - (array[i - BEAT_MS10] * 2) + array[i - 2 * BEAT_MS10]);
        sum = sum + aveBuff[avePtr];
        avePtr = avePtr + 1;
        if (avePtr == AVELENGTH)
            avePtr = 1;
        //end
        if ((i < (FIDMARK - BEAT_MS50)) || (i > (FIDMARK + BEAT_MS110)))
            if (sum > maxNoiseAve)
                maxNoiseAve = sum;
        //end
        //end

        i = i + 1;
    }//end

    if ((qrsMax - qrsMin) >= 4) {
        returnValue = ((maxNoiseAve * (50.0f / AVELENGTH)) / ((qrsMax - qrsMin) / 4));
        return returnValue;
    } else {
        // returnValue=0;
        return 0;
    }// end

}



//   //////////////////////////////////////////////////////////////////////

int searchPeakRoute(const float *ECGBuffer, int searchindex) {

    int QRSDelayLength = 100;
    // maxQueueLength=4000;
    int ModNumber = maxQueueLength;
    int maxNegative = -32768;
    int moveCount = 0;
    int movePointer = 0;
    int searchLength = 120;


    int index = searchindex;
    int currentPosition = index - QRSDelayLength;
    if (currentPosition < 1)
        currentPosition = currentPosition + maxQueueLength;
    // end
    float currentMax = (float) maxNegative;//%-32768;
    // %    currentPosition=currentPosition-frontSearch;
    int startPosition = currentPosition;
    float goalMax = (ECGBuffer[currentPosition]);//%%abs(y(currentPosition));  BeatBuffer global ECGBuffer;
    float maxQRS = (ECGBuffer[currentPosition]);
    float minQRS = (ECGBuffer[currentPosition]);
    while (moveCount < searchLength) {
        if (fabsf(goalMax) > currentMax) {
            movePointer = moveCount;
            currentMax = fabsf(goalMax);
        }//end
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (ECGBuffer[currentPosition] > maxQRS)
            maxQRS = ECGBuffer[currentPosition];
        else if (ECGBuffer[currentPosition] < minQRS)
            minQRS = ECGBuffer[currentPosition];
        //end
        // end
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currentPosition = currentPosition + 1;
        if (currentPosition == ModNumber + 1)
            currentPosition = 1;
        // end
        goalMax = ECGBuffer[currentPosition];
        moveCount = moveCount + 1;
    } //end

    int RRInterval = startPosition + movePointer;
    if (RRInterval > ModNumber)
        RRInterval = RRInterval - ModNumber;
    else if (RRInterval == ModNumber + 1)
        RRInterval = 1;
    // end
    // end
    //int    peakLostPosition=RRInterval;//%startPosition+movePointer

    return RRInterval;


}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//   float utility::data[10]={0};
//   float utility::y0=0;
//   float utility::y1=0;
//   float utility::y2=0;

float IntegerLowPass(IntegerLowPass_CTX *ctx, float datum, int init) {
//float IntegerLowPass(float datum, int init) {
//
//    static float data[10];
//    static float (*ctx_y0);
//    static float y1;
//    static float (*ctx_y2);

    float *data = ctx->data;
    float *ctx_y0 = &(ctx->y0);
    float *ctx_y1 = &(ctx->y1);
    float *ctx_y2 = &(ctx->y2);

    if (init) {
        int ii = 0;
        while (ii < 10) {
            data[ii] = 0;
            ii++;
        }
        (*ctx_y0) = 0;
        (*ctx_y1) = 0;
        (*ctx_y2) = 0;
        return 0;
    }

    int mpa = 5;
    (*ctx_y0) = 2 * (*ctx_y1) - (*ctx_y2) + datum - 2 * data[mpa - 1] + data[mpa * 2 - 1];
    //% (*ctx_y0)=2*y1-(*ctx_y2)+datum-2*data2(halfPtr)+data2(ptr);
    (*ctx_y2) = (*ctx_y1);
    (*ctx_y1) = (*ctx_y0);

    float output = (*ctx_y0) / (float) (mpa * mpa);
    int jj = 2 * mpa - 2;
    while (jj >= 0) {
        data[jj + 1] = data[jj];
        jj = jj - 1;
    }
    data[0] = datum;//%X(ptr);

    return output;
}

// y0=2*y1-y2+datum-2*data(c_mpa)+data(c_mpa*2);
//   % y0=2*y1-y2+datum-2*data2(halfPtr)+data2(ptr);
//   y2=y1;
//   y1=y0;
//   output=y0/(c_mpa*c_mpa);

//    jj=2*c_mpa-1;
//    while(jj>=1)
//      data(jj+1)=data(jj);
//      jj=jj-1;
//    end
//    data(1)=datum;%X(ptr);


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//   float utility::data_V1[10]={0};
//   float utility::y0_V1=0;
//   float utility::y1_V1=0;
//   float utility::y2_V1=0;

float IntegerLowPassV1(float datum, int init) {

    static float data_V1[10];
    static float y0_V1;
    static float y1_V1;
    static float y2_V1;
    if (init) {
        int ii = 0;
        while (ii < 10) {
            data_V1[ii] = 0;
            ii++;
        }
        y0_V1 = 0;
        y1_V1 = 0;
        y2_V1 = 0;
        return 0;
    }
    int mpa = 5;
    y0_V1 = 2 * y1_V1 - y2_V1 + datum - 2 * data_V1[mpa - 1] + data_V1[mpa * 2 - 1];
    //% y0=2*y1-y2+datum-2*data2(halfPtr)+data2(ptr);
    y2_V1 = y1_V1;
    y1_V1 = y0_V1;

    float output = y0_V1 / (float) (mpa * mpa);
    int jj = 2 * mpa - 2;
    while (jj >= 0) {
        data_V1[jj + 1] = data_V1[jj];
        jj = jj - 1;
    }
    data_V1[0] = datum;//%X(ptr);

    return output;


}






//float baseLineRemovalMulti( float ecgSample ,int channel, int init)
//{

// }

// function [ output ] =IntegerLowPassMulti(datum,multiLeadNumber,init)
float IntegerLowPassMulti(IntegerLowPassMulti_CTX *ctx, float datum, int multiLeadNumber, int init) {

    double *data = ctx->data; // static double data[10];

    double (*dataMulti)[10] = ctx->dataMulti; //static double dataMulti[multiLeadNumberMax][10];
    double *y0Multi = ctx->y0Multi; // static double y0Multi[multiLeadNumberMax];
    double *y1Multi = ctx->y1Multi; // static double y1Multi[multiLeadNumberMax];
    double *y2Multi = ctx->y2Multi; // static double y2Multi[multiLeadNumberMax];

    double *ctx_y0 = &(ctx->y0);
    double *ctx_y1 = &(ctx->y1);
    double *ctx_y2 = &(ctx->y2);
    int *ctx_mpa = &(ctx->mpa);
    int *ctx_count = &(ctx->count);


    int Fs = 500;
    int Fpa = 100;

    if (init == 1) {
        int ii = 0;
        while (ii < 10) {
            data[ii] = 0;
            ii++;
        }
        (*ctx_y0) = 0;
        (*ctx_y1) = 0;
        (*ctx_y2) = 0;

        (*ctx_mpa) = (int) round((double) Fs / (double) Fpa);

        (*ctx_count) = 0;
        int kk = 0;
        while (kk < multiLeadNumberMax) {
            int jj = 0;
            while (jj < 2 * (*ctx_mpa)) {
                dataMulti[kk][jj] = 0;
                jj = jj + 1;
            }//end

            y0Multi[kk] = 0;
            y1Multi[kk] = 0;
            y2Multi[kk] = 0;

            kk = kk + 1;
        }//end

        return 0;
    }




//persistent dataMulti;
//persistent y0Multi;
//persistent y1Multi;
//persistent y2Multi;

//multiLeadNumberMax=12;

//%  Linear integer lowpass filter
//persistent data;
//persistent (*ctx_y0);
//persistent (*ctx_y1);
//persistent (*ctx_y2);
//persistent (*ctx_mpa);
//persistent (*ctx_count);



/*
if(init==0)
  (*ctx_mpa)=round(Fs/Fpa);%2
 % sampleRate=500;
 % dataLength=length(X);
 % Xpc=zeros(1,(*ctx_mpa));%dataLengt
   data=zeros(1,2*ctx_mpa);
   %data2=zeros(1,2*ctx_mpa);
   (*ctx_y0)=0;
   (*ctx_y1)=0;
   (*ctx_y2)=0;  
   (*ctx_count)=0;
   kk=1;
   while(kk<=multiLeadNumberMax)
    jj=1;
     while(jj<=2*ctx_mpa)   
      dataMulti(kk,jj)=0;
      jj=jj+1;  
     end  
     
     y0Multi(kk)=0;     
     y1Multi(kk)=0;
     y2Multi(kk)=0;
     
     kk=kk+1;
   end
   
  % ptr=1;
  return;
end  
		 */
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%   halfPtr=ptr-(*ctx_mpa);
    if ((*ctx_count) < multiLeadNumberMax) {
        (*ctx_count) = (*ctx_count) + 1;
    }//end

    if ((*ctx_count) > multiLeadNumberMax) {
        int jj = 0;
        while (jj < 2 * (*ctx_mpa)) {
            data[jj] = dataMulti[multiLeadNumber][jj];
            jj = jj + 1;
        }//end
    }//end


    if ((*ctx_count) > multiLeadNumberMax) {
        (*ctx_y0) = y0Multi[multiLeadNumber];
        (*ctx_y1) = y1Multi[multiLeadNumber];
        (*ctx_y2) = y2Multi[multiLeadNumber];
    }
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%   if(halfPtr<1)
//%       halfPtr=halfPtr+(*ctx_mpa)*2;
//%   end
//  %while(ptr<=dataLength)       
    (*ctx_y0) = 2 * (*ctx_y1) - (*ctx_y2) + datum - 2 * data[(*ctx_mpa) - 1] + data[(*ctx_mpa) * 2 - 1];
    // / % (*ctx_y0)=2*ctx_y1-(*ctx_y2)+datum-2*data2(halfPtr)+data2(ptr);
    (*ctx_y2) = (*ctx_y1);
    (*ctx_y1) = (*ctx_y0);
    float output = (float) ((*ctx_y0) / (double) ((*ctx_mpa) * (*ctx_mpa)));

    int jj = 2 * (*ctx_mpa) - 2;
    while (jj >= 0) {
        data[jj + 1] = data[jj];
        jj = jj - 1;
    }//end
    data[0] = datum;//%X(ptr);
//%      data2(ptr)=datum;
//%      ptr=ptr+1;
//%      if(ptr==(2*ctx_mpa+1))
//%         ptr=1; 
//%      end
//  % end

/*		 
%   Tpa=((*ctx_mpa)-1);
%   T=Tpa+1;                       
%   j=1;
%   while(j<=(dataLength-T))
%      Xpc(j) =Xpc(j+T);
%      j=j+1;
%   end
%  
%   j=dataLength-(T-1);
%   while(j<=dataLength)       
%     Xpc(j)=0;
%     j=j+1;  
%   end
%   
% % Xpb(1:ns-(T))=Xpb(T+1:ns);
% % Xpb(ns-(T-1):ns)=zeros(T,1);
% %  
%  x=0:1/sampleRate:dataLength/sampleRate-1/sampleRate;
%  figure(08);                           
%  plot(x*500,Xpc);
%  
 */

    jj = 0;
    while (jj < 2 * (*ctx_mpa)) {
        dataMulti[multiLeadNumber][jj] = data[jj];
        jj = jj + 1;
    }//end

    //%   kk=1;
// %   while(kk<=multiLeadNumberMax)

    y0Multi[multiLeadNumber] = (*ctx_y0);
    y1Multi[multiLeadNumber] = (*ctx_y1);
    y2Multi[multiLeadNumber] = (*ctx_y2);

    return output;

}//end




//   //////////////////////////////////////////////////////////////////////
/// \brief utility::qrsBoundary
/// \param beat
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

//  void utility::qrsBoundary(float*beat ,float Crit,float mainWaveHeight,int QRSbegin ,int QRSfinish,int*QRSPeaksPositions,float*QRSPeaks,int index)
//  {

//   int   BEAT_MS8=4;
////#define   BEAT_MS10 5;
////    int   BEAT_MS16=8;
//    int   BEAT_MS20=10;
//  //  int   BEAT_MS50=25;
//    int   BEAT_MS90=45;
////    int   BEAT_MS200=100;
//  //  int   FIDMARK=200;
//    int   BEAT_MS250=125;
//   int    BEAT_MS40=20;
//#define    arrayLimit 300
//   int    BEAT_MS500=250;
//    int   BEAT_MS6=3;

//    float  sum=0;
// //   int   Rposition=0;

//     int  j=2;


//      float firstDiff[arrayLimit];
//      float secondDiff[arrayLimit];
//      float  fnew[arrayLimit];
//      float  data[arrayLimit];
//      float  result[arrayLimit];
//      float  y[arrayLimit];


//    //  QRSPeaks=zeros(1,BEAT_MS20);
//    //  QRSPeaksPositions=zeros(1,BEAT_MS20);

//     // data=zeros(1,arrayLimit);
//     // result=zeros(1,arrayLimit);

//      int kk=0;
//      while(kk<arrayLimit)
//      {
//         firstDiff[kk] =0;
//         secondDiff[kk]=0;
//         fnew[kk]=0;
//         result[kk]=0;
//         data[kk]=0;
//         y[kk]=0;

//         kk++;
//      }

//      kk=0;
//      while(kk<2*BEAT_MS500)
//      {
//        y[kk] =beat[kk];
//        kk++;
//      }

//     // firstDiff=zeros(1,arrayLimit);
//     // secondDiff=zeros(1,arrayLimit);
////    int   SpeakPosition=0;


//       j=2;
//      while(j<=(arrayLimit-1))
//      {
//        firstDiff[j]= (beat[j+1]-beat[j-1])/2;
//        j=j+1;
//      }//end

// //     fnew=zeros(1,arrayLimit);
//      j=3;
//      while((j>2)&&(j<(arrayLimit-2)))
//      {
//        secondDiff[j]= ((-1)*beat[j-2]-2*beat[j-1]+2*beat[j+1]+beat[j+2])/8;
//        j=j+1;
//      }//end
//      j=1;
//      while(j<=arrayLimit)
//      {
//        fnew[j]=firstDiff[j]*firstDiff[j]+secondDiff[j]*secondDiff[j];
//        j=j+1;
//       }//end

//      // data=zeros(1,arrayLimit);
//      // result=zeros(1,arrayLimit);

//       int ptr=1;
//   //    float peakPosition=0;
//       float qrsMaxValue=0;
//       j=1;
//      while(j<=arrayLimit)
//       {
//         sum=sum+fnew[j];
//         sum=sum-data[ptr];
//         data[ptr]=fnew[j];
//         ptr=ptr+1;
//        if(ptr==BEAT_MS250+1)
//           ptr=1;
//        //end
//        result[j]=sum;

//         if(((sum/BEAT_MS250)+1)>qrsMaxValue)
//           qrsMaxValue=sum/BEAT_MS250;
//         //end

//        j=j+1;

//        }//end  while(j<=arrayLimit)

//      //  y=result/BEAT_MS250;

//        j=0;
//        while(j<arrayLimit)
//        {
//           y[j]=result[j]/BEAT_MS250;
//           j++;
//        }


//       int jj=FIDMARK;
//        while((jj>2)&&(fabs(y[jj]-y[jj-2])<Crit/2)) //%% for flat top
//        {
//            jj=jj-1;
//           if(fabs(y[jj])<Crit/10)
//             break;
//           //end
//        }//end
//        while((jj>6)&&((y[jj]-y[jj-6])>(BEAT_MS250/60))&&(y[jj]>(BEAT_MS250/60)))
//        {
//          jj=jj-1;
//        }//end



//        QRSbegin=jj;
//      //  QRSPeaks=zeros(1,BEAT_MS20);
//      //  QRSPeaksPositions=zeros(1,BEAT_MS20);

//         kk=0;
//        while(kk<BEAT_MS20)
//        {
//           QRSPeaks[kk]=0;
//           QRSPeaksPositions[kk]=0;
//          kk++;
//        }

//        index=0;
//        jj=FIDMARK;
//        while(((y[jj])<(qrsMaxValue-1))&&(jj<arrayLimit))
//          jj=jj+1;
//         //end
//         QRSfinish=jj;
//         float temp;
//      //    y=beat;

//       int   i=(FIDMARK-BEAT_MS90);
//       int   peakIndex=0;

////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// while ((i<(BEAT_MS500+15))&&((i+BEAT_MS10)<arrayLimit))
//      {
//        if((i<(FIDMARK+5))&&(i>8))
//            {
//              if((fabs(y[i]-y[i-BEAT_MS10])>Crit/2)&&(fabs(y[i]-y[i+BEAT_MS10])>Crit/2))
//                  {
//                   if(((((y[i]-y[i-BEAT_MS10]))>0)&&((y[i]-y[i+BEAT_MS10])>0)) ||(((y[i]-y[i-BEAT_MS10])<0)&&(y[i]-y[i+BEAT_MS10]<0)))

//                    {
//                    index=index+1;

//                     //  downward
//                  kk=i;
//                if(y[kk]<(y[kk-BEAT_MS10]+y[kk+BEAT_MS10])/2)
//                {
//                   jj=kk-BEAT_MS10;
//                 float  temp=y[jj];
//                   peakIndex=jj;
//                   jj=jj+1;
//                   while(jj<(kk+BEAT_MS10))
//                   {
//                      if(temp>y[jj])
//                      {
//                        temp=y[jj];
//                        peakIndex=jj;
//                      }
//                    // end
//                    jj=jj+1;
//                   }//end
//                }
//               else                          //   upward
//               {
//                if((y[kk]>(y[kk-BEAT_MS10]+y[kk+BEAT_MS10])/2)&&(y[kk]>y[kk-2*BEAT_MS10]))
//                  {
//                   jj=kk-BEAT_MS10;
//                   temp=y[jj];
//                   peakIndex=jj;
//                   jj=jj+1;
//                   while(jj<(kk+BEAT_MS10))
//                   {
//                     if(temp<y[jj])
//                     {
//                       temp=y[jj];
//                       peakIndex=jj;
//                     }//end
//                    jj=jj+1;
//                   }// end
//                }// end


//              } // end

//         // %%%%%%%%%%%%%%%%%% mainWaveHeight
//             if((peakIndex>=1)&&(fabs(y[peakIndex])>mainWaveHeight*1/12))//% 22:38 28/07 2018
//                {
//                  QRSPeaksPositions[index]=peakIndex ;
//                  QRSPeaks[index]=temp;
//                }
//              else
//              index=index-1;
//                  // end

//                      i=i+2*BEAT_MS10;
//              }// end
//             }//end   if((abs(y[i]-y[i-BEAT_MS10])>C
//           } // end if((i<(FIDMARK+5))&&(i>8))
//    else //%%%%%%%%%
//               //  %%%% i>FIDMARK
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        {
//          if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
//            {
//              if((fabs(y[i]-y[i-BEAT_MS8])>Crit/2)&&(fabs(y[i]-y[i+BEAT_MS8])>Crit/2))
//               {
//                if(((((y[i]-y[i-BEAT_MS8]))>0)&&((y[i]-y[i+BEAT_MS8])>0)) ||(((y[i]-y[i-BEAT_MS8])<0)&&(y[i]-y[i+BEAT_MS8]<0)))

//                  {
//                    index=index+1;
//                   // downward
//                    kk=i;
//               if(y[kk]<(y[kk-BEAT_MS8]+y[kk+BEAT_MS8])/2)
//               {
//                     jj=kk-BEAT_MS8;
//                     temp=y[jj];
//                     peakIndex=jj;
//                     jj=jj+1;
//                     while(jj<(kk+BEAT_MS8))
//                     {
//                      if(temp>y[jj])
//                      {
//                        temp=y[jj];
//                        peakIndex=jj;
//                      }
//                      //end
//                      jj=jj+1;
//                     }// end
//               }
//              else
//                   //   upward
//               {
//                 if((y[kk]>(y[kk-BEAT_MS8]+y[kk+BEAT_MS8])/2)&&(y[kk]>y[kk+(int)(1.5*BEAT_MS20)])) // %% avoid little ripple
//                    jj=kk-BEAT_MS8;
//                    temp=y[jj];
//                    peakIndex=jj;
//                    jj=jj+1;
//                   while(jj<(kk+BEAT_MS8))
//                   {
//                     if(temp<y[jj])
//                      {
//                       temp=y[jj];
//                       peakIndex=jj;
//                      }
//                     //end
//                    jj=jj+1;
//                   }
//                   //end
//                //end

//            }//  end
//      //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//              // int ll=peakIndex;
//                 if((index>1)&&(fabs(y[peakIndex])>=mainWaveHeight*1/3)&&(peakIndex!=(QRSPeaksPositions[index-1])))//%0.58)
//                 {
//                    QRSPeaksPositions[index]=peakIndex ;
//                    QRSPeaks[index]=temp;
//                 }
//                 else
//                     index=index-1;
//                 //end
//    //  %               i=peakIndex+2*BEAT_MS6;
//                   i=i+2*BEAT_MS6;
//               } //end
//           } // end

//       else //%%%%%%%% end for  if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
//   //    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//          if((i>BEAT_MS10)&&(fabs(y[i]-y[i-BEAT_MS10])>Crit/2)&&(fabs(y[i]-y[i+BEAT_MS10])>Crit/2))
//              {
//              if((((((y[i]-y[i-BEAT_MS10]))>0)&&((y[i]-y[i+BEAT_MS10]))>0))||(((y[i]-y[i-BEAT_MS10])<0)&&(y[i]-y[i+BEAT_MS10]<0)))

//                {
//                 index=index+1;
//                  kk=i;
//                //  downward
//                if((y[kk]<(y[kk-2*BEAT_MS10]+y[kk+2*BEAT_MS10])/2)&&(y[kk]<y[kk-2*BEAT_MS20]))
//                   {
//                    jj=kk-BEAT_MS10;
//                   temp=y[jj];
//                   peakIndex=jj;
//                   jj=jj+1;
//                   while(jj<(kk+BEAT_MS10))
//                     {
//                      if(temp>y[jj])
//                       {
//                        temp=y[jj];
//                        peakIndex=jj;
//                       } //end
//                      jj=jj+1;
//                    } // end
//                 }
//             else
//              //   upward
//               {
//                 if((y[kk]>(y[kk-2*BEAT_MS10]+y[kk+2*BEAT_MS10])/2)&&(y[kk]>y[kk+2*BEAT_MS20]))//  %% avoid little ripple
//                  {
//                   jj=kk-BEAT_MS10;
//                   temp=y[jj];
//                   peakIndex=jj;
//                   jj=jj+1;
//                   while(jj<(kk+BEAT_MS10))
//                   {
//                      if(temp<y[jj])
//                      {
//                       temp=y[jj];
//                       peakIndex=jj;
//                      }//end
//                     jj=jj+1;
//                   }//end
//                 }//end
//              }//end

//            //   %  if(abs(y(peakIndex))>abs(y(FIDMARK)/2))
//                if((index>1)&&(fabs(y[peakIndex])>=mainWaveHeight*1/4)&&(peakIndex!=QRSPeaksPositions[index-1]))//%0.58)  22:38  28/07 2018
//                {
//                  QRSPeaksPositions[index]=peakIndex ;
//                  QRSPeaks[index]=temp;
//                }
//                else
//                {
//                  index=index-1;
//                }

//                i=i+2*BEAT_MS10;
//               }// end
//             }//end  line 629    if((i>BEAT_MS10)&&(fabs(y[i]-y[i-BEAT_MS10])>Crit/2)&&(fabs(y[i]-y[i+BEAT_MS10])>Crit/2))

//           }//end  line 561

//        } //end    line 563  if((i>4)&&(i>=(FIDMARK-BEAT_MS40))&&(i<=(FIDMARK+BEAT_MS40)))
//             i=i+1;




//      } // end  while ((i<(BEAT_MS500+15))&&((i+BEAT_MS10)<arrayLimit))


//  }

////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end  for void utility::qrsBoundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
