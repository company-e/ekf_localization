#ifndef INC_EKF_LOC
#define INC_EKF_LOC

#define N (4)

typedef float matNN[N][N];
typedef float vecN[N];

extern const matNN EYE;
//= {{1.0, 0.0, 0.0, 0.0},
//                   {0.0, 1.0, 0.0, 0.0},
//                   {0.0, 0.0, 1.0, 0.0},
//                   {0.0, 0.0, 0.0, 1.0}};
//
extern void ekf_estimation( vecN  xEst,
        matNN PEst,
        vecN  zObserved,
        vecN  u,
        float dt );
#endif
