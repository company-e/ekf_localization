#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ekf_localization.h"


void main()
{
    vecN  xEst = {0.0};
    matNN PEst = {0.0};
    vecN  zObserved = {0.0};
    vecN  u = {1.0, 0.1, 0.0,0.0};
    float dt = 0.10;

    memcpy(PEst, EYE, sizeof(EYE));

    for( int step = 0; step < 50/dt; step++) {
        printf("step:%4d ", step);
        for(int i = 0; i<N; i++) {
            printf("%6.3f ", xEst[i]);
        }
        printf("\n");

        zObserved[0] += u[0] * dt * cosf(xEst[2]+u[1]*dt);
        zObserved[1] += u[0] * dt * sinf(xEst[2]+u[1]*dt);//0;//u[0] * dt;

        ekf_estimation( xEst, PEst, zObserved, u, dt );


    }



}
