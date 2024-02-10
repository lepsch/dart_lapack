import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dpoequb(N, final Matrix<double> A, final int LDA, S, SCOND, AMAX, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, N;
      double             AMAX, SCOND;
      double             A( LDA, * ), S( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double             SMIN, BASE, TMP;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT, LOG, INT

      // Test the input parameters.

      // Positive definite only performs 1 pass of equilibration.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('DPOEQUB', -INFO );
         return;
      }

      // Quick return if possible.

      if ( N == 0 ) {
         SCOND = ONE;
         AMAX = ZERO;
         return;
      }

      BASE = dlamch( 'B' );
      TMP = -0.5 / LOG ( BASE );

      // Find the minimum and maximum diagonal elements.

      S[1] = A( 1, 1 );
      SMIN = S( 1 );
      AMAX = S( 1 );
      for (I = 2; I <= N; I++) { // 10
         S[I] = A( I, I );
         SMIN = min( SMIN, S( I ) );
         AMAX = max( AMAX, S( I ) );
      } // 10

      if ( SMIN <= ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 20
            if ( S( I ) <= ZERO ) {
               INFO = I;
               return;
            }
         } // 20
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 30
            S[I] = BASE ** INT( TMP * LOG( S( I ) ) );
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I)).

         SCOND = sqrt( SMIN ) / sqrt( AMAX );
      }

      }
