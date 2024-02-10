import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtrttp(final int UPLO, final int N, final Matrix<double> A, final int LDA, final int AP, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N, LDA;
      double             A( LDA, * ), AP( * );
      // ..

      bool               LOWER;
      int                I, J, K;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA

      // Test the input parameters.

      INFO = 0;
      LOWER = lsame( UPLO, 'L' );
      if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DTRTTP', -INFO );
         return;
      }

      if ( LOWER ) {
         K = 0;
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               K = K + 1;
               AP[K] = A( I, J );
            }
         }
      } else {
         K = 0;
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               K = K + 1;
               AP[K] = A( I, J );
            }
         }
      }


      }
