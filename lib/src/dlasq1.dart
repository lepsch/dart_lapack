import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlasq1(N, D, E, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N;
      double             D( * ), E( * ), WORK( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                I, IINFO;
      double             EPS, SCALE, SAFMIN, SIGMN, SIGMX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('DLASQ1', -INFO );
         return;
      } else if ( N == 0 ) {
         return;
      } else if ( N == 1 ) {
         D[1] = ( D( 1 ) ).abs();
         return;
      } else if ( N == 2 ) {
         dlas2(D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX );
         D[1] = SIGMX;
         D[2] = SIGMN;
         return;
      }

      // Estimate the largest singular value.

      SIGMX = ZERO;
      for (I = 1; I <= N - 1; I++) { // 10
         D[I] = ( D( I ) ).abs();
         SIGMX = max( SIGMX, ( E( I ) ).abs() );
      } // 10
      D[N] = ( D( N ) ).abs();

      // Early return if SIGMX is zero (matrix is already diagonal).

      if ( SIGMX == ZERO ) {
         dlasrt('D', N, D, IINFO );
         return;
      }

      for (I = 1; I <= N; I++) { // 20
         SIGMX = max( SIGMX, D( I ) );
      } // 20

      // Copy D and E into WORK (in the Z format) and scale (squaring the
      // input data makes scaling by a power of the radix pointless).

      EPS = dlamch( 'Precision' );
      SAFMIN = dlamch( 'Safe minimum' );
      SCALE = sqrt( EPS / SAFMIN );
      dcopy(N, D, 1, WORK( 1 ), 2 );
      dcopy(N-1, E, 1, WORK( 2 ), 2 );
      dlascl('G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO );

      // Compute the q's and e's.

      for (I = 1; I <= 2*N - 1; I++) { // 30
         WORK[I] = WORK( I )**2;
      } // 30
      WORK[2*N] = ZERO;

      dlasq2(N, WORK, INFO );

      if ( INFO == 0 ) {
         for (I = 1; I <= N; I++) { // 40
            D[I] = sqrt( WORK( I ) );
         } // 40
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
      } else if ( INFO == 2 ) {

      // Maximum number of iterations exceeded.  Move data from WORK
      // into D and E so the calling subroutine can try to finish

         for (I = 1; I <= N; I++) {
            D[I] = sqrt( WORK( 2*I-1 ) );
            E[I] = sqrt( WORK( 2*I ) );
         }
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO );
      }

      }
