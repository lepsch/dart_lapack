import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlagtf(N, A, LAMBDA, B, C, TOL, D, IN, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             LAMBDA, TOL;
      // ..
      // .. Array Arguments ..
      int                IN( * );
      double             A( * ), B( * ), C( * ), D( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                K;
      double             EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
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
      if ( N < 0 ) {
         INFO = -1;
         xerbla('DLAGTF', -INFO );
         return;
      }

      if (N == 0) return;

      A[1] = A( 1 ) - LAMBDA;
      IN[N] = 0;
      if ( N == 1 ) {
         if( A( 1 ) == ZERO ) IN( 1 ) = 1;
         return;
      }

      EPS = dlamch( 'Epsilon' );

      TL = max( TOL, EPS );
      SCALE1 = ( A( 1 ) ).abs() + ( B( 1 ) ).abs();
      for (K = 1; K <= N - 1; K++) { // 10
         A[K+1] = A( K+1 ) - LAMBDA;
         SCALE2 = ( C( K ) ).abs() + ( A( K+1 ) ).abs();
         if( K < ( N-1 ) ) SCALE2 = SCALE2 + ( B( K+1 ) ).abs();
         if ( A( K ) == ZERO ) {
            PIV1 = ZERO;
         } else {
            PIV1 = ( A( K ) ).abs() / SCALE1;
         }
         if ( C( K ) == ZERO ) {
            IN[K] = 0;
            PIV2 = ZERO;
            SCALE1 = SCALE2;
            if[K < ( N-1 ) ) D( K] = ZERO;
         } else {
            PIV2 = ( C( K ) ).abs() / SCALE2;
            if ( PIV2 <= PIV1 ) {
               IN[K] = 0;
               SCALE1 = SCALE2;
               C[K] = C( K ) / A( K );
               A[K+1] = A( K+1 ) - C( K )*B( K );
               if[K < ( N-1 ) ) D( K] = ZERO;
            } else {
               IN[K] = 1;
               MULT = A( K ) / C( K );
               A[K] = C( K );
               TEMP = A( K+1 );
               A[K+1] = B( K ) - MULT*TEMP;
               if ( K < ( N-1 ) ) {
                  D[K] = B( K+1 );
                  B[K+1] = -MULT*D( K );
               }
               B[K] = TEMP;
               C[K] = MULT;
            }
         }
         if( ( max( PIV1, PIV2 ) <= TL ) && ( IN( N ) == 0 ) ) IN( N ) = K;
      } // 10
      if( ( ( A( N ) ).abs() <= SCALE1*TL ) && ( IN( N ) == 0 ) ) IN( N ) = N;

      return;
      }
