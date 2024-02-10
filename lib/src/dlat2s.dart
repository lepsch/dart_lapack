import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlat2s(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Matrix<double> SA, final int LDSA, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDSA, N;
      double               SA( LDSA, * );
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- bool               lsame;
      // EXTERNAL SLAMCH, lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      RMAX = SLAMCH( 'O' );
      UPPER = lsame( UPLO, 'U' );
      if ( UPPER ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( ( A( I, J ) < -RMAX ) || ( A( I, J ) > RMAX ) ) {
                  INFO = 1;
                  GO TO 50;
               }
               SA[I][J] = double( A( I, J ) );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               if ( ( A( I, J ) < -RMAX ) || ( A( I, J ) > RMAX ) ) {
                  INFO = 1;
                  GO TO 50;
               }
               SA[I][J] = double( A( I, J ) );
            } // 30
         } // 40
      }
      } // 50

      }
