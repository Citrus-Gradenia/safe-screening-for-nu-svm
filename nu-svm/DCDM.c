/* 
 * File: DCDM.c 
 *  
 * MATLAB Coder version            : 4.1 
 * C/C++ source code generated on  : 10-May-2023 13:43:02 
 */

/* Include Files */ 
#include <math.h>
#include "DCDM.h"
#include "DCDM_emxutil.h"
#include "rand.h"
#include "DCDM_data.h"

/* Function Definitions */ 

/* 
 * 求解 1/2*alpha'*H*alpha-e'*alpha 
 *  e 是全 1 的向量 
 *  lb,ub 变量的第二个范围 
 *  n 是 alpha 元素的个数 
 *  eps 是一个极小的数 
 *  iter 是最大迭代次数 
 * Arguments    : const emxArray_real_T *H 
 *                const double lb_data[] 
 *                const int lb_size[1] 
 *                const double ub_data[] 
 *                const int ub_size[1] 
 *                double nu 
 *                double n 
 *                double eps 
 *                double iter 
 *                double alpha_data[] 
 *                int alpha_size[1] 
 * Return Type  : void 
 */
void DCDM(const emxArray_real_T *H, const double lb_data[], const int lb_size[1], const double ub_data[], const int ub_size[1], double nu, double n, double eps, double iter, double alpha_data[], int alpha_size[1])
{
    unsigned int r;
    int mti;
    emxArray_real_T *r0;
    double G;
    int i0;
    emxArray_real_T *a;
    int k;
    int i;
    boolean_T exitg1;
    int loop_ub;
    int i1;
    double b_a;
    double d0;
    (void)lb_size;
    (void)ub_size;
    r = 2U;
    state[0] = 2U;
    for (mti = 0; mti < 623; mti++) {
        r = ((r ^ r >> 30U) * 1812433253U + mti) + 1U;
        state[mti + 1] = r;
    }
    emxInit_real_T(&r0, 1);
    state[624] = 624U;
    G = nu / n;
    b_rand(n, r0);
    mti = (int)n;
    alpha_size[0] = mti;
    for (i0 = 0; i0 < mti; i0++) {
        alpha_data[i0] = G + (ub_data[i0] - G) * r0->data[i0];
    }
    emxFree_real_T(&r0);
    /*  产生范围在(nu/n*I,ub)的随机数 */
    /*  停止准则：前后两次的alpha相差很小 或 达到最大迭代次数 */
    i0 = (int)iter;
    emxInit_real_T(&a, 2);
    for (k = 0; k < i0; k++) {
        /* count<n && k<iter */
        i = 0;
        exitg1 = false;
        while ((!exitg1) && (i <= mti - 1)) {
            loop_ub = H->size[1];
            i1 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = loop_ub;
            emxEnsureCapacity_real_T(a, i1);
            for (i1 = 0; i1 < loop_ub; i1++) {
                a->data[i1] = H->data[i + H->size[0] * i1];
            }
            i1 = H->size[1];
            if ((i1 == 1) || (mti == 1)) {
                G = 0.0;
                loop_ub = a->size[1];
                for (i1 = 0; i1 < loop_ub; i1++) {
                    G += a->data[i1] * alpha_data[i1];
                }
            } else {
                G = 0.0;
                loop_ub = a->size[1];
                for (i1 = 0; i1 < loop_ub; i1++) {
                    G += a->data[i1] * alpha_data[i1];
                }
            }
            i1 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = mti;
            emxEnsureCapacity_real_T(a, i1);
            for (i1 = 0; i1 < mti; i1++) {
                a->data[i1] = 1.0;
            }
            if ((a->size[1] == 1) || (mti == 1)) {
                b_a = 0.0;
                loop_ub = a->size[1];
                for (i1 = 0; i1 < loop_ub; i1++) {
                    b_a += a->data[i1] * alpha_data[i1];
                }
            } else {
                b_a = 0.0;
                loop_ub = a->size[1];
                for (i1 = 0; i1 < loop_ub; i1++) {
                    b_a += a->data[i1] * alpha_data[i1];
                }
            }
            if (fabs(alpha_data[i] - fmax(lb_data[i], (nu - b_a) + alpha_data[i])) < eps) {
                G = fmin(G, 0.0);
            }
            /*  */
            if (fabs(alpha_data[i] - ub_data[i]) < eps) {
                G = fmax(G, 0.0);
            }
            /* 计算梯度 */
            d0 = fabs(G);
            if (d0 > eps) {
                i1 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = mti;
                emxEnsureCapacity_real_T(a, i1);
                for (i1 = 0; i1 < mti; i1++) {
                    a->data[i1] = 1.0;
                }
                if ((a->size[1] == 1) || (mti == 1)) {
                    b_a = 0.0;
                    loop_ub = a->size[1];
                    for (i1 = 0; i1 < loop_ub; i1++) {
                        b_a += a->data[i1] * alpha_data[i1];
                    }
                } else {
                    b_a = 0.0;
                    loop_ub = a->size[1];
                    for (i1 = 0; i1 < loop_ub; i1++) {
                        b_a += a->data[i1] * alpha_data[i1];
                    }
                }
                alpha_data[i] = fmin(fmax(alpha_data[i] - G / H->data[i + H->size[0] * i], fmax(lb_data[i], (nu - b_a) + alpha_data[i])), ub_data[i]);
            }
            /* 更新alpha(i) */
            /* fprintf('Iter %3d: %.4f\n', k, abs(G));  */
            if (d0 <= eps) {
                exitg1 = true;
            } else {
                /* 不更新aplha */
                i++;
            }
        }
    }
    emxFree_real_T(&a);
}
/* 
 * File trailer for DCDM.c 
 *  
 * [EOF] 
 */
