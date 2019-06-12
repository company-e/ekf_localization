#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "ekf_localization.h"

const matNN EYE = {{1.0, 0.0, 0.0, 0.0},
                   {0.0, 1.0, 0.0, 0.0},
                   {0.0, 0.0, 1.0, 0.0},
                   {0.0, 0.0, 0.0, 1.0}};

const matNN Q   = {{0.1 * 0.1, 0.0,       0.0,                             0.0},
                   {0.0,       0.1 * 0.1, 0.0,                             0.0},
                   {0.0,       0.0,       (1.0*M_PI/180) * (1.0*M_PI/180), 0.0},
                   {0.0,       0.0,       0.0,                             1.0}};

const matNN R   = {{1.0*1.0, 0.0,     0.0, 0.0},
                   {0.0,     1.0*1.0, 0.0, 0.0},
                   {0.0,     0.0,     0.0, 0.0},
                   {0.0,     0.0,     0.0, 0.0}};

void mat_trans( matNN b, const matNN a )
{
    int i, j;

    for ( i = 0; i<N; i++) {
        for ( j = 0; j<N; j++ ) {
            b[i][j] = a[j][i];
        }
    }
}

void multply_mat_mat( matNN c, const matNN a, const matNN b )
{
    int i, j, k;
    float s;

    for ( i = 0; i < N; i++ ) {
        for ( j = 0; j < N; j++ ) {
            s = 0;
            for ( k = 0; k < N; k++) {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
}

void multply_mat_vec( vecN c, const matNN a, const vecN b )
{
    for( int i = 0; i<N; i++ ) {
        c[i] = 0;
        for( int j = 0; j<N; j++ ) {
            c[i] += a[i][j] * b[j];
        }
    }
}

void multply_vec_mat( vecN c, const vecN a, const matNN b )
{
    for( int i = 0; i<N; i++ ) {
        c[i] = 0;
        for( int j = 0; j<N; j++ ) {
            c[i] += a[j] * b[j][i];
        }
    }
}

void add_mat( matNN c, const matNN a, const matNN b )
{
    int i, j;

    for ( i = 0; i < N; i++ ) {
        for ( j = 0; j < N; j++ ) {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
}

void add_vec( vecN c, const vecN a, const vecN b )
{
    int i;
    for (i = 0; i<N; i++) {
        c[i] = a[i] + b[i];
    }
}

void sub_mat( matNN c, const matNN a, const matNN b )
{
    int i, j;

    for ( i = 0; i < N; i++ ) {
        for ( j = 0; j < N; j++ ) {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
}

void sub_vec( vecN c, const vecN a, const vecN b )
{
    int i;
    for (i = 0; i<N; i++ ){
        c[i] = a[i] - b[i];
    }
}

void matrix_inv_2by2( matNN b, matNN a )
{
    float tmp;
    //matNN a_inv;
    tmp = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if( fabsf(tmp) > 0.0000001 ) {
        b[0][0] =      (1/tmp) * a[1][1];
        b[0][1] = -1 * (1/tmp) * a[0][1];
        b[1][0] = -1 * (1/tmp) * a[1][0];
        b[1][1] =      (1/tmp) * a[0][0];
    } else {
        b = NULL;
    }
}

void jacobH( matNN jH )
{
    jH[0][0] = 1.0;
    jH[0][1] = 0.0;
    jH[0][2] = 0.0;
    jH[0][3] = 0.0;

    jH[1][0] = 0.0;
    jH[1][1] = 1.0;
    jH[1][2] = 0.0;
    jH[1][3] = 0.0;

}

void jacobH_T( matNN jH_T )
{
    jH_T[0][0] = 1.0;
    jH_T[1][0] = 0.0;
    jH_T[2][0] = 0.0;
    jH_T[3][0] = 0.0;

    jH_T[0][1] = 0.0;
    jH_T[1][1] = 1.0;
    jH_T[2][1] = 0.0;
    jH_T[3][1] = 0.0;

}

void motion_model( vecN xpred, const vecN x,  const vecN u, const float dt )
{
    /* ( x, y ,yaw ) */
    /* x[0] : px */
    /* x[1] : py */
    /* x[2] : yaw */
    /* x[3] : v */
    /* u[0] : v */
    /* u[1] : yawrate */
    /* xpred = F@x +B@u*/
    matNN F = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0}};

    float cos_yaw = cosf(x[2]);
    float sin_yaw = sinf(x[2]);

    matNN B = {
        {dt*cos_yaw, 0.0, 0.0, 0.0},
        {dt*sin_yaw, 0.0, 0.0, 0.0},
        {0.0,        dt,  0.0, 0.0},
        {1.0,        0.0, 0.0, 0.0}};
    vecN F_x = {0.0};
    vecN B_u = {0.0};

    multply_mat_vec( F_x, F, x);
    multply_mat_vec( B_u, B, u);
    add_vec( xpred, F_x, B_u ); 

}

void jacobF(matNN jF, float dt, vecN x, vecN u)
{
    /* 
    Jacobian of Motion Model

    motion model
    x_{t+1} = x_t+v*dt*cos(yaw)
    y_{t+1} = y_t+v*dt*sin(yaw)
    yaw_{t+1} = yaw_t+omega*dt
    v_{t+1} = v{t}
    so
    dx/dyaw = -v*dt*sin(yaw)
    dx/dv = dt*cos(yaw)
    dy/dyaw = v*dt*cos(yaw)
    dy/dv = dt*sin(yaw)
    */
    float yaw = x[2];
    float vlc = u[0];
    float sin_yaw = sinf(yaw);
    float cos_yaw = cosf(yaw);
    dt = 0.1;
    memcpy(jF, EYE, sizeof(EYE));
    jF[0][2] = -dt * vlc * sin_yaw;
    jF[1][2] =  dt * vlc * cos_yaw;
    jF[0][3] =  dt * cos_yaw;
    jF[1][3] =  dt * sin_yaw;
    jF[0][0] = 1.0;
    jF[1][1] = 1.0;
    jF[2][2] = 1.0;
    jF[3][3] = 1.0;
}

void PredictState( vecN xpred, const vecN x, const vecN u, const float dt )
{
   motion_model( xpred, x, u, dt ); 
}

void PredictErrCovarince( matNN PPred, const matNN jF, const matNN PEst )
{
    /* PPred = jF@PEst@jF.T + Q */
    matNN jFT         = {0.0};
    matNN jF_PEst     = {0.0};
    matNN jF_PEst_jFT = {0.0};

    mat_trans( jFT, jF );
    multply_mat_mat( jF_PEst, jF, PEst );
    multply_mat_mat( jF_PEst_jFT, jF_PEst, jFT );
    add_mat( PPred, jF_PEst_jFT, Q );

}

void CalcKalmanGain( matNN K, matNN PPred, matNN jH )
{
    matNN jH_PPred     = {0.0};
    matNN jH_PPred_jHT = {0.0};
    matNN jH_T         = {0.0};
    matNN S            = {0.0};
    matNN S_inv        = {0.0};
    matNN PPred_jHT    = {0.0};

    jacobH_T( jH_T );

    /* S = jH@PPred@jH.T + R */
    multply_mat_mat( jH_PPred, jH, PPred);
    multply_mat_mat( jH_PPred_jHT, jH_PPred, jH_T );
    add_mat( S, jH_PPred_jHT, R);

    /* K = PPred@jH.T@np.linalg.inv(S) */
    matrix_inv_2by2( S_inv, S );
    multply_mat_mat( PPred_jHT, PPred, jH_T );
    multply_mat_mat( K, PPred_jHT, S_inv );

}

void UpdateState(
        vecN xEst,
        vecN xPred,
        matNN K,
        vecN zObserved )
{
    vecN zPred = {0.0};
    vecN  y     = {0.0};
    vecN  K_y   = {0.0};

    /* xEst = xPred + K@y */
    zPred[0] = xPred[0];
    zPred[1] = xPred[1];

    y[0] = zObserved[0] - zPred[0];
    y[1] = zObserved[1] - zPred[1];

    multply_mat_vec( K_y, K, y );
    add_vec( xEst, xPred, K_y );

}

void UpdateErrCovariance(
        matNN PEst,
        matNN K,
        matNN jH,
        matNN PPred )
{
    matNN K_jH       = {0.0};
    matNN eye_K_jH   = {0.0};


    /* PEst = (np.eye(len(xEst)) - K@jH)@PPred */
    multply_mat_mat( K_jH, K, jH);
    sub_mat( eye_K_jH, EYE, K_jH);
    multply_mat_mat( PEst, eye_K_jH, PPred );
}

/*
 * xEst[0] : px
 * xEst[1] : py
 * xEst[2] : yaw
 * xEst[3] : v
 * zObserved[0] : px
 * zObserved[1] : py
 * u[0]         : vlc
 * u[1]         : yawrate
 *
 */
void ekf_estimation(
        vecN  xEst,
        matNN PEst,
        vecN  zObserved,
        vecN  u,
        float dt )
{
    vecN xPred  = {0.0, 0.0,0.0,0.0}; /* 4行1列 */
    matNN PPred = {0.0}; /* 4行4列 */
    matNN jF    = {0.0}; /* 4行4列 */
    matNN jH    = {0.0}; /* 2行4列 */
    matNN K     = {0.0}; /* 4行2列*/

    /* Predict Step */
    PredictState( xPred, xEst, u, dt );
    jacobF( jF, dt, xPred, u );
    PredictErrCovarince( PPred, jF, PEst );

    /* Update Step */
    jacobH( jH );
    CalcKalmanGain( K, PPred, jH );
    UpdateState( xEst, xPred, K, zObserved );
    UpdateErrCovariance( PEst, K, jH, PPred );
}
