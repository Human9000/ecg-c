//
// Created by 陈昞翱 on 2023/5/13.
//

#ifndef ECG_ALGORITHM_C_UTIL_H
#define ECG_ALGORITHM_C_UTIL_H

#include <stdio.h>

void m_log_int(int data) {
    FILE *out = NULL;
    if (!(out = fopen("../log/log.txt", "a"))) {
        return;
    }

    fprintf(out, "%d\n", data);
    fclose(out);
}

void m_log_float(float data) {
    FILE *out = NULL;
    if (!(out = fopen("../log/log.txt", "a"))) {
        return;
    }

    fprintf(out, "%f\n", data);
    fclose(out);
}

#endif //ECG_ALGORITHM_C_UTIL_H
