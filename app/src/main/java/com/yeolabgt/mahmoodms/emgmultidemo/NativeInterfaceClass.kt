package com.yeolabgt.mahmoodms.emgmultidemo

/**
 * Created by Musa Mahmood on 1/14/2018.
 *
 */

class NativeInterfaceClass {


    @Throws(IllegalArgumentException::class)
    external fun jmainInitialization(b: Boolean): Int
    @Throws(IllegalArgumentException::class)
    external fun jTrainingRoutineKNN2(data: DoubleArray): DoubleArray
    @Throws(IllegalArgumentException::class)
    external fun jClassifyUsingKNNv4(data: DoubleArray, kparams: DoubleArray, k: Double): Double

    companion object {
        init {
            System.loadLibrary("ssvep-lib")
        }
    }
}