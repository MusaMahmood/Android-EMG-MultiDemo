//
// Created by mahmoodms on 4/3/2017.
//

#include "rt_nonfinite.h"
#include "ctrainingRoutineKNN2.h"
#include "classifyArmEMG4.h"

/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define  LOG_TAG "jniExecutor-cpp"
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

extern "C" {
JNIEXPORT jdouble JNICALL
/**
 *
 * @param env
 * @param jobject1
 * @param allData 750x3 vector of data
 * @param params KNN features 1x9920
 * @param LastY
 * @return
 */
Java_com_yeolabgt_mahmoodms_emgmultidemo_NativeInterfaceClass_jClassifyUsingKNNv4(
        JNIEnv *env, jobject jobject1, jdoubleArray allData, jdoubleArray params, jdouble knn) {
    jdouble *X1 = env->GetDoubleArrayElements(allData, nullptr);
    jdouble *PARAMS = env->GetDoubleArrayElements(params, nullptr);
    if (X1 == nullptr) LOGE("ERROR - C_ARRAY IS NULL");
    return classifyArmEMG4(X1, PARAMS, knn);
}
}

extern "C" {
JNIEXPORT jdoubleArray JNICALL
Java_com_yeolabgt_mahmoodms_emgmultidemo_NativeInterfaceClass_jTrainingRoutineKNN2(JNIEnv *env, jobject obj, jdoubleArray allData) {
    jdouble *X = env->GetDoubleArrayElements(allData, nullptr);
    if (X== nullptr) LOGE("ERROR - C_ARRAY");
    double KNN_PARAMS[9920];
    //TODO: INSERT ANALYSIS FUNCTION HERE
    ctrainingRoutineKNN2(X,KNN_PARAMS);
    jdoubleArray mReturnArray = env->NewDoubleArray(9920);
    env->SetDoubleArrayRegion(mReturnArray, 0, 9920, &KNN_PARAMS[0]);
    return mReturnArray;
}
}

extern "C" {
JNIEXPORT jint JNICALL
Java_com_yeolabgt_mahmoodms_emgmultidemo_NativeInterfaceClass_jmainInitialization(
        JNIEnv *env, jobject obj, jboolean terminate) {
    if (!(bool) terminate) {
        ctrainingRoutineKNN2_initialize();
        return 0;
    } else {
        return -1;
    }
}
}
