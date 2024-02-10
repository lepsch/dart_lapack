import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgttrf(N, DL, D, DU, DU2, IPIV, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N;
      int                IPIV( * );
      double             D( * ), DL( * ), DU( * ), DU2( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                I;
      double             FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('DGTTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize IPIV(i) = i and DU2(I) = 0

      for (I = 1; I <= N; I++) { // 10
         IPIV[I] = I;
      } // 10
      for (I = 1; I <= N - 2; I++) { // 20
         DU2[I] = ZERO;
      } // 20

      for (I = 1; I <= N - 2; I++) { // 30
         if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {

            // No row interchange required, eliminate DL(I)

            if ( D( I ) != ZERO ) {
               FACT = DL( I ) / D( I );
               DL[I] = FACT;
               D[I+1] = D( I+1 ) - FACT*DU( I );
            }
         } else {

            // Interchange rows I and I+1, eliminate DL(I)

            FACT = D( I ) / DL( I );
            D[I] = DL( I );
            DL[I] = FACT;
            TEMP = DU( I );
            DU[I] = D( I+1 );
            D[I+1] = TEMP - FACT*D( I+1 );
            DU2[I] = DU( I+1 );
            DU[I+1] = -FACT*DU( I+1 );
            IPIV[I] = I + 1;
         }
      } // 30
      if ( N > 1 ) {
         I = N - 1;
         if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {
            if ( D( I ) != ZERO ) {
               FACT = DL( I ) / D( I );
               DL[I] = FACT;
               D[I+1] = D( I+1 ) - FACT*DU( I );
            }
         } else {
            FACT = D( I ) / DL( I );
            D[I] = DL( I );
            DL[I] = FACT;
            TEMP = DU( I );
            DU[I] = D( I+1 );
            D[I+1] = TEMP - FACT*D( I+1 );
            IPIV[I] = I + 1;
         }
      }

      // Check for a zero on the diagonal of U.

      for (I = 1; I <= N; I++) { // 40
         if ( D( I ) == ZERO ) {
            INFO = I;
            GO TO 50;
         }
      } // 40
      } // 50

      }
