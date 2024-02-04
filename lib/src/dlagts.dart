import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlagts(JOB, N, A, B, C, D, IN, Y, TOL, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, JOB, N;
      double             TOL;
      // ..
      // .. Array Arguments ..
      int                IN( * );
      double             A( * ), B( * ), C( * ), D( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                K;
      double             ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( ( ( JOB ).abs() > 2 ) || ( JOB == 0 ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('DLAGTS', -INFO );
         return;
      }

      if (N == 0) return;

      EPS = dlamch( 'Epsilon' );
      SFMIN = dlamch( 'Safe minimum' );
      BIGNUM = ONE / SFMIN;

      if ( JOB < 0 ) {
         if ( TOL <= ZERO ) {
            TOL = ( A( 1 ) ).abs();
            if (N > 1) TOL = max( TOL, ( A( 2 ) ).abs(), ( B( 1 ) ) ).abs();
            for (K = 3; K <= N; K++) { // 10
               TOL = max( TOL, ( A( K ) ).abs(), ( B( K-1 ) ).abs(), ( D( K-2 ) ) ).abs();
            } // 10
            TOL = TOL*EPS;
            if (TOL == ZERO) TOL = EPS;
         }
      }

      if ( ( JOB ).abs() == 1 ) {
         for (K = 2; K <= N; K++) { // 20
            if ( IN( K-1 ) == 0 ) {
               Y[K] = Y( K ) - C( K-1 )*Y( K-1 );
            } else {
               TEMP = Y( K-1 );
               Y[K-1] = Y( K );
               Y[K] = TEMP - C( K-1 )*Y( K );
            }
         } // 20
         if ( JOB == 1 ) {
            for (K = N; K >= 1; K--) { // 30
               if ( K <= N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 );
               } else if ( K == N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 );
               } else {
                  TEMP = Y( K );
               }
               AK = A( K );
               ABSAK = ( AK ).abs();
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ( TEMP ).abs()*SFMIN > ABSAK ) {
                        INFO = K;
                        return;
                     } else {
                        TEMP = TEMP*BIGNUM;
                        AK = AK*BIGNUM;
                     }
                  } else if ( ( TEMP ).abs() > ABSAK*BIGNUM ) {
                     INFO = K;
                     return;
                  }
               }
               Y[K] = TEMP / AK;
            } // 30
         } else {
            for (K = N; K >= 1; K--) { // 50
               if ( K <= N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 );
               } else if ( K == N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 );
               } else {
                  TEMP = Y( K );
               }
               AK = A( K );
               PERT = sign( TOL, AK );
               } // 40
               ABSAK = ( AK ).abs();
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ( TEMP ).abs()*SFMIN > ABSAK ) {
                        AK = AK + PERT;
                        PERT = 2*PERT;
                        GO TO 40;
                     } else {
                        TEMP = TEMP*BIGNUM;
                        AK = AK*BIGNUM;
                     }
                  } else if ( ( TEMP ).abs() > ABSAK*BIGNUM ) {
                     AK = AK + PERT;
                     PERT = 2*PERT;
                     GO TO 40;
                  }
               }
               Y[K] = TEMP / AK;
            } // 50
         }
      } else {

         // Come to here if  JOB = 2 or -2

         if ( JOB == 2 ) {
            for (K = 1; K <= N; K++) { // 60
               if ( K >= 3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 );
               } else if ( K == 2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 );
               } else {
                  TEMP = Y( K );
               }
               AK = A( K );
               ABSAK = ( AK ).abs();
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ( TEMP ).abs()*SFMIN > ABSAK ) {
                        INFO = K;
                        return;
                     } else {
                        TEMP = TEMP*BIGNUM;
                        AK = AK*BIGNUM;
                     }
                  } else if ( ( TEMP ).abs() > ABSAK*BIGNUM ) {
                     INFO = K;
                     return;
                  }
               }
               Y[K] = TEMP / AK;
            } // 60
         } else {
            for (K = 1; K <= N; K++) { // 80
               if ( K >= 3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 );
               } else if ( K == 2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 );
               } else {
                  TEMP = Y( K );
               }
               AK = A( K );
               PERT = sign( TOL, AK );
               } // 70
               ABSAK = ( AK ).abs();
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ( TEMP ).abs()*SFMIN > ABSAK ) {
                        AK = AK + PERT;
                        PERT = 2*PERT;
                        GO TO 70;
                     } else {
                        TEMP = TEMP*BIGNUM;
                        AK = AK*BIGNUM;
                     }
                  } else if ( ( TEMP ).abs() > ABSAK*BIGNUM ) {
                     AK = AK + PERT;
                     PERT = 2*PERT;
                     GO TO 70;
                  }
               }
               Y[K] = TEMP / AK;
            } // 80
         }

         for (K = N; K >= 2; K--) { // 90
            if ( IN( K-1 ) == 0 ) {
               Y[K-1] = Y( K-1 ) - C( K-1 )*Y( K );
            } else {
               TEMP = Y( K-1 );
               Y[K-1] = Y( K );
               Y[K] = TEMP - C( K-1 )*Y( K );
            }
         } // 90
      }
      }
