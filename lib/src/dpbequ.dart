import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dpbequ(final int UPLO, final int N, final int KD, final Matrix<double> AB, final int LDAB, final int S, final int SCOND, final int AMAX, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, N;
      double             AMAX, SCOND;
      double             AB( LDAB, * ), S( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               UPPER;
      int                I, J;
      double             SMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDAB < KD+1 ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('DPBEQU', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         SCOND = ONE;
         AMAX = ZERO;
         return;
      }

      if ( UPPER ) {
         J = KD + 1;
      } else {
         J = 1;
      }

      // Initialize SMIN and AMAX.

      S[1] = AB( J, 1 );
      SMIN = S( 1 );
      AMAX = S( 1 );

      // Find the minimum and maximum diagonal elements.

      for (I = 2; I <= N; I++) { // 10
         S[I] = AB( J, I );
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
            S[I] = ONE / sqrt( S( I ) );
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = sqrt( SMIN ) / sqrt( AMAX );
      }
      }
