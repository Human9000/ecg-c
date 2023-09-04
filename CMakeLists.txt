cmake_minimum_required(VERSION 3.24)
project(ecg_algorithm_c C)

set(CMAKE_C_STANDARD 11)

include_directories(include)

link_libraries(m)

add_executable(${PROJECT_NAME}
        src/main.c
        src/algorithm/ArrhythmiaAnalysis20230306.c
        src/algorithm/beatFeaturesGet.c
        src/algorithm/ECGPreProcess.c
        src/algorithm/pWaveDetected.c
        src/algorithm/qrsDetect20230211.c
        src/algorithm/tWaveFeaturesGet.c
        src/algorithm/utility.c
        src/algorithm/common_util.c
        )

