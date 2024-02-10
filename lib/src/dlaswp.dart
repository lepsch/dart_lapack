import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaswp(N, A, LDA, K1, K2, IPIV, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, K1, K2, LDA, N;
      int                IPIV( * );
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, I1, I2, INC, IP, IX, IX0, J, K, N32;
      double             TEMP;

      // Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
      // K1 through K2.

      if ( INCX > 0 ) {
         IX0 = K1;
         I1 = K1;
         I2 = K2;
         INC = 1;
      } else if ( INCX < 0 ) {
         IX0 = K1 + ( K1-K2 )*INCX;
         I1 = K2;
         I2 = K1;
         INC = -1;
      } else {
         return;
      }

      N32 = ( N / 32 )*32;
      if ( N32 != 0 ) {
         for (J = 1; J <= N32; J += 32) { // 30
            IX = IX0;
            for (I = I1; INC < 0 ? I >= I2 : I <= I2; I += INC) { // 20
               IP = IPIV( IX );
               if ( IP != I ) {
                  for (K = J; K <= J + 31; K++) { // 10
                     TEMP = A( I, K );
                     A[I][K] = A( IP, K );
                     A[IP][K] = TEMP;
                  } // 10
               }
               IX = IX + INCX;
            } // 20
         } // 30
      }
      if ( N32 != N ) {
         N32 = N32 + 1;
         IX = IX0;
         for (I = I1; INC < 0 ? I >= I2 : I <= I2; I += INC) { // 50
            IP = IPIV( IX );
            if ( IP != I ) {
               for (K = N32; K <= N; K++) { // 40
                  TEMP = A( I, K );
                  A[I][K] = A( IP, K );
                  A[IP][K] = TEMP;
               } // 40
            }
            IX = IX + INCX;
         } // 50
      }

      }
