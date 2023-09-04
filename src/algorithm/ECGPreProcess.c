#include <math.h>
#include "ECGPreProcess.h"

double a[7] = {1, -5.11025492374927, 10.9229366097765, -12.5120053127799, 8.10246853688899, -2.80937380187589,
               0.406229011001932};
double b[7] = {0.070577093298849442, -0.24761750225933435, 0.28350475666195968, -0.0000000000000003134252559660019,
               -0.28350475666195912, 0.24761750225933404, -0.070577093298849314};


float IIRBandFilter(IIRBandFilter_CTX *ctx, float inputdata, int initial) {
//float IIRBandFilter( float inputdata, int initial) {
//    volatile static float IIRpy[iirFilterLength];
//    volatile static float IIRpx[iirFilterLength];

    float *IIRpy = ctx->IIRpy;
    float *IIRpx = ctx->IIRpx;


    if (initial == 1) {
        for (int ii = 0; ii < iirFilterLength; ii++) {
            IIRpy[ii] = 0;
            IIRpx[ii] = 0;
        }

        return 0;
    }

    IIRpx[0] = inputdata;
    volatile float sum = 0;

//          for jj=1:1:7
//            sum=sum+b(jj)*IIRpx(jj);
//          end
    int jj = 0;
    while (jj < iirFilterLength) {
        sum = sum + (float) b[jj] * IIRpx[jj];
        jj = jj + 1;
    }

//          for jj=2:1:7
//            sum=sum-a(jj)*IIRpy(jj);
//          end
    jj = 1;
    while (jj < iirFilterLength) {
        sum = sum - (float) a[jj] * IIRpy[jj];
        jj = jj + 1;
    }



//          k=7;
//          while (k>=2)
//              IIRpx(k)=IIRpx(k-1);
//              k=k-1;
//          end


    int k = iirFilterLength - 1;
    while (k >= 1) {
        IIRpx[k] = IIRpx[k - 1];
        k = k - 1;
    }

//           k=7;
//          while (k>=3)
//             IIRpy(k)=IIRpy(k-1);
//              k=k-1;
//          end

    k = iirFilterLength - 1;
    while (k >= 2) {
        IIRpy[k] = IIRpy[k - 1];
        k = k - 1;
    }



//        IIRpy(2)=sum;
    IIRpy[1] = sum;
//       % yy(ii)=sum;
//       % sum=0.0;
//       %xx(ii)=yy(ii);
//      % y(ii)=yy(ii);
//        outdata=sum;
    //        outputData=sum;


    return sum;


}

//

void dataSort(const float *inputArray, float *outputArray, int length) {

    int index_i = 0;
    while (index_i < length) {
        outputArray[index_i] = inputArray[index_i];
        index_i++;
    }

    index_i = 0;
    while (index_i < (length - 1)) {
        int index_j = index_i + 1;
        while (index_j < (length)) {
            if (outputArray[index_j] < outputArray[index_i]) {
                float temp = outputArray[index_i];
                outputArray[index_i] = outputArray[index_j];
                outputArray[index_j] = temp;
            }

            index_j = index_j + 1;
        }
        index_i = index_i + 1;
    }


}


void binaryInsertSort(float *dataArray, float currentData, float data, int lenght) {


    if (currentData == data)
        return;

    int mid;
    int begin = 0;
    int finish = lenght - 1;
    while (begin <= finish) {
        mid = ((begin + finish) / 2);
        if (currentData == dataArray[mid])
            break;
        else if (dataArray[mid] < currentData)
            begin = mid + 1;
        else if (dataArray[mid] > currentData)
            finish = mid - 1;

    }

//    %% From the deleted position to move data array to left
    int j = mid;
    while (j < lenght - 1) {
        dataArray[j] = dataArray[j + 1];
        j = j + 1;
    }

//    %%%%%%%%%%%%%%%%%%%%%
//    %deletePoint=mid;
//    %%%%%%%%%%%%%%%%%%%%
//    % %% Search the deleted data position
    begin = 0;
    finish = lenght - 1;
//    %finish=length;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (begin <= finish) {
        mid = ((begin + finish) / 2);
        if (dataArray[mid] < data)
            begin = mid + 1;
        else
            finish = mid - 1;

    }

//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    %insertPoint=finish+1;
//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    % %% From the deleted position to move data array to left

    j = lenght - 2;
    while (j > (finish)) {
        dataArray[j + 1] = dataArray[j];
        j = j - 1;
    }

    if ((finish + 1) <= (lenght - 1))
        dataArray[finish + 1] = data;
    else
        dataArray[finish] = data;


}


float baseLineRemoval(baseLineRemoval_CTX *ctx, float inputdata, int initial) {
    float *resultOne = ctx->resultOne;
    float *templateOne = ctx->templateOne;
    float *resultTwo = ctx->resultTwo;
    float *templateTwo = ctx->templateTwo;
    float *realTimeData = ctx->realTimeData;
    int *ctx_count = &(ctx->count);


    if (initial == 1) {

        // initialize the static array
        (*ctx_count) = 0;
        int jj = 0;
        while (jj < lenghtLevelOne) {
            resultOne[jj] = 0;
            templateOne[jj] = 0;
            jj++;
        }


        jj = 0;
        while (jj < lengthLevelTwo) {
            resultTwo[jj] = 0;
            templateTwo[jj] = 0;
            jj++;
        }

        jj = 0;
        while (jj < (lengthLevelTwo - 1) / 2) {
//        while (jj < (lengthLevelTwo)) {
            realTimeData[jj] = 0;
            jj++;
        }

        return 0;
    }


    (*ctx_count) = (*ctx_count) + 1;
    if ((*ctx_count) > (lengthLevelTwo + lenghtLevelOne)) {
        (*ctx_count) = (lengthLevelTwo + lenghtLevelOne) + 10;
    }

    float tempBuffer;
    float outputdata = 0;
    int index_A = 0;
    tempBuffer = templateOne[index_A];
    float baseWanderData;

    while (index_A < (lenghtLevelOne - 1)) {
        templateOne[index_A] = templateOne[index_A + 1];
        index_A = index_A + 1;
    }

    templateOne[lenghtLevelOne - 1] = inputdata;
    if ((*ctx_count) == lenghtLevelOne) {
        dataSort(templateOne, resultOne, lenghtLevelOne);
        baseWanderData = resultOne[(lenghtLevelOne - 1) / 2];
    } else if ((*ctx_count) > lenghtLevelOne) {

        binaryInsertSort(resultOne, tempBuffer, inputdata, lenghtLevelOne);
        baseWanderData = resultOne[(lenghtLevelOne - 1) / 2 - 1];

    }


    if ((*ctx_count) >= lenghtLevelOne) {
        int index_B = 0;
        tempBuffer = templateTwo[1];
        while (index_B < (lengthLevelTwo - 1)) {
            templateTwo[index_B] = templateTwo[index_B + 1];
            index_B = index_B + 1;
        }
        templateTwo[lengthLevelTwo - 1] = baseWanderData;


        if ((*ctx_count) == (lengthLevelTwo + lenghtLevelOne - 1)) {

            dataSort(templateTwo, resultTwo, lengthLevelTwo);
            baseWanderData = resultTwo[(lengthLevelTwo - 1) / 2];
            outputdata = realTimeData[1] - baseWanderData;
        } else {

            index_B = 0;
            while (index_B < ((lengthLevelTwo - 1) / 2 - 1)) {
                realTimeData[index_B] = realTimeData[index_B + 1];
                index_B = index_B + 1;
            }

            realTimeData[(lengthLevelTwo - 1) / 2 - 1] = inputdata;

            if ((*ctx_count) >= (lengthLevelTwo + lenghtLevelOne)) {
                binaryInsertSort(resultTwo, tempBuffer, baseWanderData, lengthLevelTwo);
                baseWanderData = resultTwo[(lengthLevelTwo - 1) / 2];
                outputdata = realTimeData[1] - baseWanderData;

            } else
                outputdata = 0;
        }

    }

    return outputdata;

}



//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

float baseLineRemovalMulti(baseLineRemovalMulti_CTX *ctx, float ecgSample, int multiLeadNumber, int init) {
    //int multiLeadNumberMax=12;
    //%% Multi leads
    float (*resultOneMulti)[lenghtLevelOne] = ctx->resultOneMulti; //static float resultOneMulti[multiLeadNumberMax][lenghtLevelOne];
    float (*templateOneMulti)[lenghtLevelOne] = ctx->templateOneMulti; //static float templateOneMulti[multiLeadNumberMax][lenghtLevelOne];
    float (*resultTwoMulti)[lengthLevelTwo] = ctx->resultTwoMulti; //static float resultTwoMulti[multiLeadNumberMax][lengthLevelTwo];
    float (*templateTwoMulti)[lengthLevelTwo] = ctx->templateTwoMulti; //static float templateTwoMulti[multiLeadNumberMax][lengthLevelTwo];
    float (*realTimeDataMulti)[(lengthLevelTwo - 1) /
                               2] = ctx->realTimeDataMulti; //static float realTimeDataMulti[multiLeadNumberMax][(lengthLevelTwo - 1) / 2];
    int *countMulti = ctx->countMulti; // static int countMulti[12];


    float *resultOne = ctx->resultOne; // static float resultOne[lenghtLevelOne];
    float *templateOne = ctx->templateOne; // static float templateOne[lenghtLevelOne];
    float *resultTwo = ctx->resultTwo; // static float resultTwo[lengthLevelTwo];
    float *templateTwo = ctx->templateTwo; // static float templateTwo[lengthLevelTwo];
    float *realTimeData = ctx->realTimeData; // static float realTimeData[(lengthLevelTwo - 1) / 2];

    int *ctx_count = &(ctx->count);
// resultOneMulti;
//persistent resultTwoMulti;


    float y = ecgSample;

    if (init == 1) {
        (*ctx_count) = 0;

        // initialize the static array
        (*ctx_count) = 0;
        int jj = 0;
        while (jj < lenghtLevelOne) {
            resultOne[jj] = 0;
            templateOne[jj] = 0;
            jj++;
        }


        jj = 0;
        while (jj < lengthLevelTwo) {
            resultTwo[jj] = 0;
            templateTwo[jj] = 0;
            jj++;
        }

        jj = 0;
        while (jj < (lengthLevelTwo - 1) / 2) {
//        while (jj < (lengthLevelTwo)) {
            realTimeData[jj] = 0;
            jj++;
        }


        jj = 0;
        while (jj < multiLeadNumberMax) {

            countMulti[jj] = 0;
            int kk = 0;
            while (kk < lenghtLevelOne) {
                resultOneMulti[jj][kk] = 0;
                templateOneMulti[jj][kk] = 0;
                kk = kk + 1;
            }//end

            kk = 0;
            while (kk < lengthLevelTwo) {
                resultTwoMulti[jj][kk] = 0;
                templateTwoMulti[jj][kk] = 0;
                kk = kk + 1;
            }//end

            kk = 0;
            while (kk < (lengthLevelTwo - 1) / 2) {
                realTimeDataMulti[jj][kk] = 0;
                kk = kk + 1;
            }//end

            jj = jj + 1;
        }// end

        return 0;

    }//end  if(init==0)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    if (countMulti[multiLeadNumber] < constValue) {
        countMulti[multiLeadNumber] = countMulti[multiLeadNumber] + 1;
    }//end
    (*ctx_count) = countMulti[multiLeadNumber];

    int jj = 0;
    while (((*ctx_count) > lenghtLevelOne) && (jj < lenghtLevelOne)) {
        resultOne[jj] = resultOneMulti[multiLeadNumber][jj];
        templateOne[jj] = templateOneMulti[multiLeadNumber][jj];
        jj = jj + 1;
    }//end

    jj = 0;
    while (((*ctx_count) > lenghtLevelOne) && (jj < lengthLevelTwo)) {
        resultTwo[jj] = resultTwoMulti[multiLeadNumber][jj];
        templateTwo[jj] = templateTwoMulti[multiLeadNumber][jj];
        jj = jj + 1;
    }//end

    jj = 0;
    while (((*ctx_count) > lenghtLevelOne) && (jj < (lengthLevelTwo - 1) / 2)) {
        realTimeData[jj] = realTimeDataMulti[multiLeadNumber][jj];
        jj = jj + 1;
    }//end


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
//%   index_B=1;
//%   while(index_B<(lengthLevelTwo))
//%     realTimeData(index_B)=realTimeData(index_B+1);
//%     index_B=index_B+1;
//%   end
//%     realTimeData(lengthLevelTwo)=data;

    float outputdata = 0;
    float baseWanderData;
    int index_A = 0;
    float tempBuffer = templateOne[1];
    while (index_A < lenghtLevelOne) {
        templateOne[index_A] = templateOne[index_A + 1];
        index_A = index_A + 1;
    }//end
    //% tempBuffer=templateOne(1);
    templateOne[lenghtLevelOne - 1] = ecgSample;
    if ((*ctx_count) == lenghtLevelOne) {
        //resultOne = DataSort(templateOne,lenghtLevelOne-1 );
        //baseWanderData=resultOne[(lenghtLevelOne-1)/2];
        dataSort(templateOne, resultOne, lenghtLevelOne);
        baseWanderData = resultOne[(lenghtLevelOne - 1) / 2];
    } else {
        if ((*ctx_count) > lenghtLevelOne) {
            //        resultOne  = BinaryInsertSort( resultOne,tempBuffer,ecgSample,lenghtLevelOne-1 );
            //%    resultOne  = QuickInsertSort( resultOne,tempBuffer,data,lenghtLevelOne );
            //         baseWanderData=resultOne[(lenghtLevelOne-1)/2];

            binaryInsertSort(resultOne, tempBuffer, ecgSample, lenghtLevelOne);
            baseWanderData = resultOne[(lenghtLevelOne - 1) / 2 - 1];

        }//end
    }//end

    int index_B;
    if ((*ctx_count) >= lenghtLevelOne - 1) {
        index_B = 0;
        tempBuffer = templateTwo[1];
        while (index_B < lengthLevelTwo) {
            templateTwo[index_B] = templateTwo[index_B + 1];
            index_B = index_B + 1;
        }//end
        templateTwo[lengthLevelTwo - 1] = baseWanderData;

        if ((*ctx_count) == (lengthLevelTwo + lenghtLevelOne - 1)) {
            //          resultTwo = DataSort(templateTwo,lengthLevelTwo );
            //          baseWanderData=resultTwo((lengthLevelTwo-1)/2);
            //        %  y=realTimeData((lengthLevelTwo+lenghtLevelOne-2)/2+1)-baseWanderData;
            //         y=realTimeData(1)-baseWanderData;
            dataSort(templateTwo, resultTwo, lengthLevelTwo);
            baseWanderData = resultTwo[(lengthLevelTwo - 1) / 2];
            outputdata = realTimeData[1] - baseWanderData;

        } else {
            index_B = 0;
            //%     while(index_B<(lengthLevelTwo))
            while (index_B < ((lengthLevelTwo - 1) / 2)) {
                realTimeData[index_B] = realTimeData[index_B + 1];
                //if((index_B+1)>=175)
                //   kkkk=0;
                //end
                index_B = index_B + 1;
            }//end
            //%  realTimeData(lengthLevelTwo)=data;
            realTimeData[(lengthLevelTwo - 1) / 2 - 1] = ecgSample;

            if ((*ctx_count) >= (lengthLevelTwo + lenghtLevelOne)) {
                // resultTwo  = BinaryInsertSort( resultTwo,tempBuffer,baseWanderData,lengthLevelTwo);
                //%              resultTwo  = QuickInsertSort( resultTwo,tempBuffer,baseWanderData,lengthLevelTwo);
                // baseWanderData=resultTwo((lengthLevelTwo-1)/2);
//%               y=realTimeData((lengthLevelTwo-1)/2)-baseWanderData;
                //  y=realTimeData(1)-baseWanderData;

                binaryInsertSort(resultTwo, tempBuffer, baseWanderData, lengthLevelTwo);
                baseWanderData = resultTwo[(lengthLevelTwo - 1) / 2];
                outputdata = realTimeData[1] - baseWanderData;

            }//end

        }//end

//%      if((*ctx_count)<=398)
//%         y=y-data; 
//%      end

    }//end
//%  if((*ctx_count)==401)
//%  llll=0;
//%  end
    if ((*ctx_count) < 402)
        //  %  y=y-data;
        y = 0;
    //end  



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
//   int kk=0;
//  while(kk<=multiLeadNumberMax)
    //   %countMulti(multiLeadNumber)=(*ctx_count);
//    kk=kk+1;
//  end

    jj = 0;
    while (jj < lenghtLevelOne) {
        resultOneMulti[multiLeadNumber][jj] = resultOne[jj];
        templateOneMulti[multiLeadNumber][jj] = templateOne[jj];
        jj = jj + 1;
    }//end

    int ll = 0;
    while (ll < lengthLevelTwo) {
        resultTwoMulti[multiLeadNumber][ll] = resultTwo[ll];
        templateTwoMulti[multiLeadNumber][ll] = templateTwo[ll];
        ll = ll + 1;
    }//end

    int pp = 0;
    while (pp < (lengthLevelTwo - 1) / 2) {
        realTimeDataMulti[multiLeadNumber][pp] = realTimeData[pp];
        pp = pp + 1;
    }//end

    return outputdata;
}


//  ??????????????????????????????????????????????????????????

//  50HZ notch  filter
float NotchFilter_50hz(NotchFilter_50hz_CTX *ctx,float datum, int channelNumber, int init) {
    //	int  sampleRate=500;
    //  int  notchRate=50;
    //	int  cycleLength=sampleRate/notchRate;
    //	static sampleRate;
    //  static notchRate;
    //  static cycleLength;
    float (*notchQueue)[cycleLength + 2] = ctx->notchQueue; //static float notchQueue[ChannelMax][cycleLength + 2];
    float (*notchData)[cycleLength] = ctx->notchData; //static float notchData[ChannelMax][cycleLength];
    float (*pbuffer)[cycleLength] = ctx->pbuffer; //static float pbuffer[ChannelMax][cycleLength];
    int *index = ctx->index; // static int index[ChannelMax];

    int        *ctx_flag = &(ctx->flag);
    int        *ctx_c_threshold = &(ctx->c_threshold);

    if (init == 1) {
        // sampleRate=500;
        // notchRate=50;
        (*ctx_c_threshold) = 18;
        (*ctx_flag) = 1;
        int jj = 0;
        while (jj < ChannelMax) {
            int ii = 0;
            while (ii < cycleLength) {
                notchData[jj][ii] = 0;
                pbuffer[jj][ii] = 0;
                notchQueue[jj][ii] = 0;
                ii = ii + 1;
            }

            notchQueue[jj][cycleLength] = 0;
            notchQueue[jj][cycleLength + 1] = 0;
            index[jj] = 0;
            jj = jj + 1;
        }

        // notchQueue[cycleLength] = 0;
        // notchQueue[cycleLength + 1] = 0;
        // index = 0;
        //  notchQueue=zeros(1,cycleLength+2);
        //  notchData =zeros(1,cycleLength);
        //  pbuffer=zeros(1,cycleLength);

        //%     //50HZ removal
        //%         int (*ctx_c_threshold) = 18;
        //%         int sampleRate = 500;
        //%         int notchRate = 50;
        //%         int (*ctx_flag) = 1;
        //%         int index = 1;
        //%         float[] notchQueue = new float[12];
        //%         float[] notchData = new float[10];
        //%         float[] pbufer = new float[10];
        return 0;
    }
//end
/*

% %(*ctx_flag)=1;
% %j=1;
% 
% for ii=1:dataLength%30000
%    data(ii)=0;
% end;
% for ii=1:10
%    p(ii)=0;
% end;
%   
%   
% for ii=1:1:dataLength-cycleLength-2
%    if j==cycleLength+1;
%       j=1;
%    end;
%  DifferA=y(ii+cycleLength)-y(ii);
%  DifferB=y(ii+cycleLength+1)-y(ii+1);
%  if abs(DifferA-DifferB)<(*ctx_c_threshold)
%     (*ctx_flag)=(*ctx_flag)-1;
%     if (*ctx_flag)==0
%           (*ctx_flag)=1;
%           temp=0;
%              for jj=0:1:cycleLength-1
%                temp=temp+y(ii+jj);     
%              end;
%          differ=y(cycleLength+ii)-y(ii);
%          data(ii+cycleLength/2-1)=(temp-differ/2)/cycleLength;
%          p(j)=y(ii+cycleLength/2-1)-data(ii+cycleLength/2-1);
%     else
%       data(ii+cycleLength/2-1)=y(ii+cycleLength/2-1)-p(j);
%     end;
%   else
%    (*ctx_flag)=cycleLength;
%    data(ii+cycleLength/2-1)=y(ii+cycleLength/2-1)-p(j);
% end;      
% 
% j=j+1;
% 
% end;
% 
% figure(28);
% plot(x*360,data);
% grid on;
% y=data;
% figure(66);
% fs=360;%500;%pp;%2400;
% N=10000;%nn;%6000;
% n=0:N-1;
% t=n/fs;
% z=fft(data,N);
% mag=abs(z);
% f=n*fs/N;
% plot(f(1:N/2),mag(1:N/2),'r--');
% title('Signal Amplitude');
% xlabel('Frequency(Hz)');
% title(' Spectrum Analysis After 50HZ Removal  Filter ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%         public float NotchRejection(float data)
%         {
%             //50HZ removal
%             int i;
%             float temp;
%             int cycleLength = sampleRate / notchRate;
%             float differA, differB, differ;
% 
%             // notchQueue[cycleLength] = data;
%             //  notchQueue[cycleLength+1] = data;
%             for (i = 0; i < cycleLength + 1; i++)
%                 notchQueue[i] = notchQueue[i + 1];
% 
%             if (index == cycleLength)
%                 index = 1;
%                differA = notchQueue[cycleLength - 1] - notchQueue[0];
%                differB = notchQueue[cycleLength] - notchQueue[1];
%             //   differA = notchQueue[cycleLength ] - notchQueue[0];
%             //   differB = notchQueue[cycleLength+1] - notchQueue[1];
%             if (System.Math.Abs(differA - differB) < (*ctx_c_threshold))
%             {
%                 (*ctx_flag) -= 1;
%                 if ((*ctx_flag) == 0)
%                 {
%                     (*ctx_flag) = 1;
%                     temp = 0;
%                     
%                     for (i = 0; i < cycleLength; i++)
%                     {
%                         temp += notchQueue[i];
%                     }
%                     differ = notchQueue[cycleLength - 1] - notchQueue[0];
%                     notchData[cycleLength / 2 - 1] = (temp - differ) / cycleLength;
%                     pbufer[index] = notchQueue[cycleLength / 2 - 1] - notchData[cycleLength / 2 - 1];
% 
%                 }
%                 else
%                     notchData[cycleLength / 2 - 1] = notchQueue[cycleLength / 2 - 1] - pbufer[index];
% 
% 
%             }
%             else
%                 (*ctx_flag) = cycleLength;
%             notchData[cycleLength / 2 - 1] = notchQueue[cycleLength / 2 - 1] - pbufer[index];
%             temp = notchQueue[cycleLength / 2 - 1] - pbufer[index];
%             
% 
%             index ++;
%             //50HZ notch removal
% 
%             return temp;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*/
    notchQueue[channelNumber][cycleLength + 1] = datum;
/// for i=1:1:cycleLength+1
    //  notchQueue(i) = notchQueue(i + 1);
    //end
    int ii = 0;
    while (ii <= cycleLength) {
        notchQueue[channelNumber][ii] = notchQueue[channelNumber][ii + 1];
        ii = ii + 1;
    }
    // if (index == cycleLength)
    //    index = 1;
    // end
    if (index[channelNumber] == cycleLength - 1) {
        index[channelNumber] = 0;
    }


    //    DifferA = notchQueue(cycleLength ) - notchQueue(1);
    //    DifferB = notchQueue(cycleLength+1) - notchQueue(2);

    float DifferA = notchQueue[channelNumber][cycleLength - 1] - notchQueue[channelNumber][0];
    float DifferB = notchQueue[channelNumber][cycleLength] - notchQueue[channelNumber][1];

    if (fabsf(DifferA - DifferB) < (float) (*ctx_c_threshold)) {
        (*ctx_flag) = (*ctx_flag) - 1;
        if ((*ctx_flag) == 0) {
            (*ctx_flag) = 1;
            float temp = 0;
            //   for i = 1:1: cycleLength
            //    temp = temp+notchQueue(i);
            //  end

            ii = 0;
            while (ii < cycleLength) {
                temp = temp + notchQueue[channelNumber][ii];
                ii = ii + 1;
            }
            //  differ = notchQueue(cycleLength ) - notchQueue(1);
            //  notchData(cycleLength / 2 - 1) = (temp - differ) / cycleLength;
            //  pbuffer(index) = notchQueue(cycleLength / 2 - 1) - notchData(cycleLength / 2 - 1);
            float differ = notchQueue[channelNumber][cycleLength - 1] - notchQueue[channelNumber][0];
            notchData[channelNumber][cycleLength / 2 - 1] = (temp - differ) / cycleLength;
            pbuffer[channelNumber][index[channelNumber]] =
                    notchQueue[channelNumber][cycleLength / 2 - 1] - notchData[channelNumber][cycleLength / 2 - 1];
        } else {
            notchData[channelNumber][cycleLength / 2 - 1] =
                    notchQueue[channelNumber][cycleLength / 2 - 1] - pbuffer[channelNumber][index[channelNumber]];
        }
        //  end
    } else {
        (*ctx_flag) = cycleLength - 1;
    }

    //end
    //  notchData(cycleLength / 2 - 1) = notchQueue(cycleLength / 2 - 1) - pbuffer(index);
    //  temp = notchQueue(cycleLength / 2 - 1) - pbuffer(index);
    //  index =index+1;
    //  y= temp;

    notchData[channelNumber][cycleLength / 2 - 1] =
            notchQueue[channelNumber][cycleLength / 2 - 1] - pbuffer[channelNumber][index[channelNumber]];
    float temp = notchQueue[channelNumber][cycleLength / 2 - 1] - pbuffer[channelNumber][index[channelNumber]];
    index[channelNumber] = index[channelNumber] + 1;

    //   y= temp;

    return temp;


}





// //
//  float  IntegerLowPass(float inputdata,int initial){
//
//		static float datap[sampleRate/c_fpa*2];
//    static float y0;
//    static float y1;
//    static float y2;
//		if(initial==1)
//		{
//			y0=0;
//			y1=0;
//			y2=0;
//			int ii=0;
//			while(ii<sampleRate/c_fpa*2){
//
//				datap[ii]=0;
//				ii=ii+1;
//
//			}
//
//			return 0;
//		}
//
//			 y0=2*y1-y2+inputdata-2*datap[c_mpa-1]+datap[c_mpa*2-1];
//    //   % y0=2*y1-y2+inputdata-2*data2(halfPtr)+data2(ptr);
//       y2=y1;
//       y1=y0;
//       float  output=y0/(c_mpa*c_mpa);

//       int jj=2*c_mpa-2;
//        while(jj>=0){
//          datap[jj+1]=datap[jj];
//          jj=jj-1;
//        }
//        datap[0]=inputdata;

//       return output;
//
//
//
//	}
//
//
//
//
//
//
//
//
//
//
//
 
 

