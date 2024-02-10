import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgetc2(N, A, LDA, IPIV, JPIV, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, N;
      int                IPIV( * ), JPIV( * );
      double             A( LDA, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, IP, IPV, J, JP, JPV;
      double             BIGNUM, EPS, SMIN, SMLNUM, XMAX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGER, DSWAP
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;

      // Set constants to control overflow

      EPS = dlamch( 'P' );
      SMLNUM = dlamch( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Handle the case N=1 by itself

      if ( N == 1 ) {
         IPIV[1] = 1;
         JPIV[1] = 1;
         if ( ( A( 1, 1 ) ).abs() < SMLNUM ) {
            INFO = 1;
            A[1][1] = SMLNUM;
         }
         return;
      }

      // Factorize A using complete pivoting.
      // Set pivots less than SMIN to SMIN.

      for (I = 1; I <= N - 1; I++) { // 40

         // Find max element in matrix A

         XMAX = ZERO;
         for (IP = I; IP <= N; IP++) { // 20
            for (JP = I; JP <= N; JP++) { // 10
               if ( ( A( IP, JP ) ).abs() >= XMAX ) {
                  XMAX = ( A( IP, JP ) ).abs();
                  IPV = IP;
                  JPV = JP;
               }
            } // 10
         } // 20
         if (I == 1) SMIN = max( EPS*XMAX, SMLNUM );

         // Swap rows

         if (IPV != I) dswap( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA );
         IPIV[I] = IPV;

         // Swap columns

         if (JPV != I) dswap( N, A( 1, JPV ), 1, A( 1, I ), 1 );
         JPIV[I] = JPV;

         // Check for singularity

         if ( ( A( I, I ) ).abs() < SMIN ) {
            INFO = I;
            A[I][I] = SMIN;
         }
         for (J = I + 1; J <= N; J++) { // 30
            A[J][I] = A( J, I ) / A( I, I );
         } // 30
         dger(N-I, N-I, -ONE, A( I+1, I ), 1, A( I, I+1 ), LDA, A( I+1, I+1 ), LDA );
      } // 40

      if ( ( A( N, N ) ).abs() < SMIN ) {
         INFO = N;
         A[N][N] = SMIN;
      }

      // Set last pivots to N

      IPIV[N] = N;
      JPIV[N] = N;

      }
