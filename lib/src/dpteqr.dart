import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dpteqr(COMPZ, N, D, E, Z, LDZ, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             COMPZ;
      int                INFO, LDZ, N;
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDSQR, DLASET, DPTTRF, XERBLA
      double             C( 1, 1 ), VT( 1, 1 );
      int                I, ICOMPZ, NRU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters.

      INFO = 0;

      if ( lsame( COMPZ, 'N' ) ) {
         ICOMPZ = 0;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ICOMPZ = 2;
      } else {
         ICOMPZ = -1;
      }
      if ( ICOMPZ < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ( LDZ < 1 ) || ( ICOMPZ > 0 && LDZ < max( 1, N ) ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DPTEQR', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if (ICOMPZ > 0) Z( 1, 1 ) = ONE;
         return;
      }
      if (ICOMPZ == 2) dlaset( 'Full', N, N, ZERO, ONE, Z, LDZ );

      // Call DPTTRF to factor the matrix.

      dpttrf(N, D, E, INFO );
      if (INFO != 0) return;
      for (I = 1; I <= N; I++) { // 10
         D[I] = sqrt( D( I ) );
      } // 10
      for (I = 1; I <= N - 1; I++) { // 20
         E[I] = E( I )*D( I );
      } // 20

      // Call DBDSQR to compute the singular values/vectors of the
      // bidiagonal factor.

      if ( ICOMPZ > 0 ) {
         NRU = N;
      } else {
         NRU = 0;
      }
      dbdsqr('Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO );

      // Square the singular values.

      if ( INFO == 0 ) {
         for (I = 1; I <= N; I++) { // 30
            D[I] = D( I )*D( I );
         } // 30
      } else {
         INFO = N + INFO;
      }

      return;
      }
