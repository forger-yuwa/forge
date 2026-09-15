#pragma once
// リミッタ関数 (Venkatakrishnan / Barth-Jespersen) の device 実体。limiter_d.cu からの**テキスト移動** (式・型・順序不変)。
// 受動種の無次元化リミッタ (passiveLimiter_d.cuh) と単体試験が同じ関数を使うためヘッダに置く。
// __global__ カーネルを持つ限り 1 つのライブラリ TU (limiter_d.cu) からだけ include すること (単体試験 TU は別実行体)。
#include "flowFormat.hpp"

// Limiters for Unstructured Higher-Order Accurate Solutions of the Euler Equations
// Krzysztof Michalak

__device__ flow_float venkata_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                      flow_float delta_m, flow_float volume) {

    flow_float K = 1.f;
    flow_float eps2 = K*K*K*volume;
    //return (x*x + 2.0*x + eps*eps)/(x*x + x + 2.0 + eps*eps);
    flow_float res;

    // K11: 元は /(...)/delta_m と除算2回。/(denom*delta_m) に統合して除算1回に（compute律速の limiter 向け）。
    if (delta_m > 1e-20f) {
        flow_float delta_p = delta_p_max;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0f*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else if (delta_m < -1e-20f) {
        flow_float delta_p = delta_p_min;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0f*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else {
        res = 1.0f;
    }

    return res;
}

__device__ flow_float barth_Jespersen_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                              flow_float delta_m, flow_float volume) {

    flow_float res;

    if (delta_m > 1e-20f) {
        res = min(1.0f, delta_p_max/delta_m);
    } else if (delta_m < -1e-20f) {
        res = min(1.0f, delta_p_min/delta_m);
    } else {
        res = 1.0f;
    }

    return min(res, 1.0f);
}


