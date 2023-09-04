//
// Created by 陈昞翱 on 2023/5/6.
//
#include "common_util.h"

int counterAdjust(int counter, int modeValue) {
    int tempCount = counter;
    if (tempCount > modeValue) {
        tempCount = tempCount - modeValue;
    } else if (tempCount < 0) {
        tempCount = modeValue + tempCount;
    }//end
    //end

    return tempCount;
}




