import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaqsb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, UPLO;
      int                KD, LDAB, N;
      double             AMAX, SCOND;
      double             AB( LDAB, * ), S( * );
      // ..

      double             ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      int                I, J;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Quick return if possible

      if ( N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = dlamch( 'Safe minimum' ) / dlamch( 'Precision' );
      LARGE = ONE / SMALL;

      if ( SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No equilibration

         EQUED = 'N';
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( lsame( UPLO, 'U' ) ) {

            // Upper triangle of A is stored in band format.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = max( 1, J-KD ); I <= J; I++) { // 10
                  AB[KD+1+I-J][J] = CJ*S( I )*AB( KD+1+I-J, J );
               } // 10
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               for (I = J; I <= min( N, J+KD ); I++) { // 30
                  AB[1+I-J][J] = CJ*S( I )*AB( 1+I-J, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      }
