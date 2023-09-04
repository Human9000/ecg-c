#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include "ecgmain.h"

enum Lead {
    I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
};

void get_data(char *str, float *rtn_data_arr, const int rtn_data_arr_len, const enum Lead *lead_num_arr,
              const int lead_num_arr_len) {
    assert(rtn_data_arr && rtn_data_arr_len > 0);
    assert(lead_num_arr_len >= 2);
    float lead_I = 0;
    float lead_II = 0;
    float lead_III = 0;
    float lead_aVR = 0;
    float lead_aVL = 0;
    float lead_aVF = 0;
    float lead_V1 = 0;
    float lead_V2 = 0;
    float lead_V3 = 0;
    float lead_V4 = 0;
    float lead_V5 = 0;
    float lead_V6 = 0;
    char *p;
    char *ptr;
    float data[] = {0, 0};
    for (int i = 0; i < 2; ++i) {
        assert(i < rtn_data_arr_len);
        char *data_str = strtok_r(i == 0 ? str : NULL, " \n", &p);
        if (data_str == NULL) break;
        data[i] = strtof(data_str, &ptr);
        switch (lead_num_arr[i]) {
            case I:
                lead_I = data[i];
                break;
            case II:
                lead_II = data[i];
                break;
            case III:
                lead_III = data[i];
                break;
            case aVR:
                lead_aVR = data[i];
                break;
            case aVL:
                lead_aVL = data[i];
                break;
            case aVF:
                lead_aVF = data[i];
                break;
            case V1:
                lead_V1 = data[i];
                break;
            case V2:
                lead_V2 = data[i];
                break;
            case V3:
                lead_V3 = data[i];
                break;
            case V4:
                lead_V4 = data[i];
                break;
            case V5:
                lead_V5 = data[i];
                break;
            case V6:
                lead_V6 = data[i];
                break;
        }
    }

    if (lead_I == 0) {
        lead_I = data[0];
    }
    if (lead_II == 0) {
        lead_II = data[1];
    }
    if (lead_III == 0) {
        lead_III = data[0];
    }
    if (lead_aVR == 0) {
        lead_aVR = data[1];
    }
    if (lead_aVL == 0) {
        lead_aVL = data[0];
    }
    if (lead_aVF == 0) {
        lead_aVF = data[1];
    }
    if (lead_V1 == 0) {
        lead_V1 = data[0];
    }
    if (lead_V2 == 0) {
        lead_V2 = data[1];
    }
    if (lead_V3 == 0) {
        lead_V3 = data[0];
    }
    if (lead_V4 == 0) {
        lead_V4 = data[1];
    }
    if (lead_V5 == 0) {
        lead_V5 = data[0];
    }
    if (lead_V6 == 0) {
        lead_V6 = data[1];
    }


    rtn_data_arr[0] = lead_I;
    rtn_data_arr[1] = lead_II;
    rtn_data_arr[2] = lead_III;
    rtn_data_arr[3] = lead_aVR;
    rtn_data_arr[4] = lead_aVL;
    rtn_data_arr[5] = lead_aVF;
    rtn_data_arr[6] = lead_V1;
    rtn_data_arr[7] = lead_V2;
    rtn_data_arr[8] = lead_V3;
    rtn_data_arr[9] = lead_V4;
    rtn_data_arr[10] = lead_V5;
    rtn_data_arr[11] = lead_V6;
}

void print_arrhythmia_type(const arrhythmia_type *input) {
    printf("{\n"
           "  supra_arrhythmia: {\n"
           "    apcTotal: %d,\n"
           "    japcTotal: %d,\n"
           "    doubleApc: %d,\n"
           "    apcBigeminyCount: %d,\n"
           "    apcTrigeminyCount: %d,\n"
           "    atrialTachycardiaConut: %d,\n"
           "    longestApcTachycardia: %d,\n"
           "    fastestApcTachy: %d,\n"
           "    atrial_fibrillation: %d,\n"
           "    atrial_flutter: %d\n"
           "  },\n"
           "  ventricular_arrhythmia: {\n"
           "    pvcCount: %d,\n"
           "    doublePvc: %d,\n"
           "    vpcBigeminyCount: %d,\n"
           "    vpcTrigeminyCount: %d,\n"
           "    ventricularTachycardiaCount: %d,\n"
           "    longestPvcTachycardia: %d,\n"
           "    fastestPvcTachy: %d,\n"
           "    lbbbCount: %d,\n"
           "    rbbbCount: %d\n"
           "  },\n"
           "  sinus_arrhythmia: {\n"
           "    sinusArrhythmiaType: %d,\n"
           "    sinusIrregularityType: %d,\n"
           "    sinusArrestType: %d,\n"
           "    sinusBradycardiaType: %d,\n"
           "    sinusTachycardiaType: %d,\n"
           "    sickSinusSyndromeType: %d,\n"
           "  },\n"
           "  escape_Count: {\n"
           "    junctional_escape: %d,\n"
           "    ventricular_escape: %d\n"
           "  },\n"
           "  block_Count: {\n"
           "    one_degreeCount: %d,\n"
           "    two_degreeCount: %d,\n"
           "    three_degreeCount: %d,\n"
           "    high_degreeCount: %d,\n"
           "    AsystolyCount: %d\n"
           "  },\n"
           "  st_value: {\n"
           "    st_value: %f,\n"
           "    type: %d,\n"
           "  }\n"
           "}\n",
           input->supra_arrhythmia.apcTotal,
           input->supra_arrhythmia.japcTotal,
           input->supra_arrhythmia.doubleApc,
           input->supra_arrhythmia.apcBigeminyCount,
           input->supra_arrhythmia.apcTrigeminyCount,
           input->supra_arrhythmia.atrialTachycardiaConut,
           input->supra_arrhythmia.longestApcTachycardia,
           input->supra_arrhythmia.fastestApcTachy,
           input->supra_arrhythmia.atrial_fibrillation,
           input->supra_arrhythmia.atrial_flutter,
           input->ventricular_arrhythmia.pvcCount,
           input->ventricular_arrhythmia.doublePvc,
           input->ventricular_arrhythmia.vpcBigeminyCount,
           input->ventricular_arrhythmia.vpcTrigeminyCount,
           input->ventricular_arrhythmia.ventricularTachycardiaCount,
           input->ventricular_arrhythmia.longestPvcTachycardia,
           input->ventricular_arrhythmia.fastestPvcTachy,
           input->ventricular_arrhythmia.lbbbCount,
           input->ventricular_arrhythmia.rbbbCount,
           input->sinus_arrhythmia.sinusArrhythmiaType,
           input->sinus_arrhythmia.sinusIrregularityType,
           input->sinus_arrhythmia.sinusArrestType,
           input->sinus_arrhythmia.sinusBradycardiaType,
           input->sinus_arrhythmia.sinusTachycardiaType,
           input->sinus_arrhythmia.sickSinusSyndromeType,
           input->escape_Count.junctional_escape,
           input->escape_Count.ventricular_escape,
           input->block_Count.one_degreeCount,
           input->block_Count.two_degreeCount,
           input->block_Count.three_degreeCount,
           input->block_Count.high_degreeCount,
           input->block_Count.AsystolyCount,
           input->st_value.st_value,
           input->st_value.type
    );

}

ArrhythmiaAnalysis20230311_CTX *new_ecg_ctx( ) {
    ArrhythmiaAnalysis20230311_CTX *ctx = (ArrhythmiaAnalysis20230311_CTX *) malloc(
            sizeof(ArrhythmiaAnalysis20230311_CTX));

    const float data_arr[12] = {0};
    ctx->count = 0;
    ctx->rrCount = 0;

    ArrhythmiaAnalysis20230311(ctx,
                               data_arr,
                               12,
                               &(ctx->rrCount),
                               &(ctx->ret),
                               true,
                               0);
    return ctx;
}

int free_ect_ctx(ArrhythmiaAnalysis20230311_CTX *ctx) {
    free(ctx);
    return true;
}

int ecg_run(ArrhythmiaAnalysis20230311_CTX *ctx,
            float *data_arr) {
    ctx->count = ctx->count + 1;
    return ArrhythmiaAnalysis20230311(ctx,
                                      data_arr,
                                      12,
                                      &(ctx->rrCount),
                                      &(ctx->ret),
                                      false,
                                      ctx->count);
}


arrhythmia_type ecg_get_ret(ArrhythmiaAnalysis20230311_CTX *ctx) {
    return ctx->ret;
}
int ecg_get_count(ArrhythmiaAnalysis20230311_CTX *ctx) {
    return ctx->count;
}


void run(char *input_filename, const enum Lead *lead_num_arr, const int lead_num_arr_len) {
    FILE *input = NULL;
    if (!(input = fopen(input_filename, "r"))) {
        return;
    }

    char data_str[1000];
    const int data_arr_len = 12;
    float data_arr[12];
    memset(data_arr, 0, sizeof data_arr);
    // ============== ecg use case start ===========================
    ArrhythmiaAnalysis20230311_CTX *ctx = new_ecg_ctx();  // 新建ecg需要的相关资源
    while (fgets(data_str, 1000, input)) {
        get_data(data_str, data_arr, data_arr_len, lead_num_arr, lead_num_arr_len);
        ecg_run(ctx, data_arr);// 流式输入一个长度为12的数组
    }
    arrhythmia_type ret = ecg_get_ret(ctx); // 获得返回值
    int count = ecg_get_count(ctx); // 获得资源数量
    free_ect_ctx(ctx); // 释放资源
    // ============== ecg use case end ===========================
    printf("total_data_count: %d\n", count);
    printf("%s:\n", input_filename);
    print_arrhythmia_type(&ret);
    printf("-----------------------------------\n");
    fclose(input);
}

void run_back(char *input_filename, const enum Lead *lead_num_arr, const int lead_num_arr_len) {
    FILE *input = NULL;
    if (!(input = fopen(input_filename, "r"))) {
        return;
    }

    char data_str[1000];
    arrhythmia_type ret;
    const int data_arr_len = 12;
    float data_arr[12];
    memset(data_arr, 0, sizeof data_arr);
    int rrCount = 0;
    ArrhythmiaAnalysis20230311_CTX ctx;
    ArrhythmiaAnalysis20230311(&ctx, data_arr, data_arr_len, &rrCount, &ret, true, 0);
    int count = 0;
    while (fgets(data_str, 1000, input)) {
        count++;
        get_data(data_str, data_arr, data_arr_len, lead_num_arr, lead_num_arr_len);
        ArrhythmiaAnalysis20230311(&ctx, data_arr, data_arr_len, &rrCount, &ret, false, count);
    }
    printf("total_data_count: %d\n", count);
    printf("%s:\n", input_filename);
    print_arrhythmia_type(&ret);
    printf("-----------------------------------\n");
    fclose(input);
}

int main() {
    system("chcp 65001 > nul"); // windows 下加入这行代码可以防止控制台打印中文乱码
    printf("开始执行...\n");
    clock_t start = clock();
    enum Lead lead_num_arr[] = {II, V1};
//    run_back("../data/500/103m.mat.txt", lead_num_arr, 2);
    run("../data/500/103m.mat.txt", lead_num_arr, 2);
    clock_t end = clock();
    unsigned long duration = (end - start) / CLOCKS_PER_SEC;
    printf("完成，共执行了 %lu 分 %lu 秒", duration / 60, duration % 60);
    return 0;
}


