import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlarzt(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      double             T( LDT, * ), TAU( * ), V( LDV, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                I, INFO, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DTRMV, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // Check for currently supported options

      INFO = 0;
      if ( !lsame( DIRECT, 'B' ) ) {
         INFO = -1;
      } else if ( !lsame( STOREV, 'R' ) ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('DLARZT', -INFO );
         return;
      }

      for (I = K; I >= 1; I--) { // 20
         if ( TAU( I ) == ZERO ) {

            // H(i)  =  I

            for (J = I; J <= K; J++) { // 10
               T[J][I] = ZERO;
            } // 10
         } else {

            // general case

            if ( I < K ) {

               // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**T

               dgemv('No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 );

               // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

               dtrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
            }
            T[I][I] = TAU( I );
         }
      } // 20
      }
