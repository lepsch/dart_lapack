// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cgtsv(final int N, final int NRHS, final int DL, final int D, final int DU, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final B = B_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDB, N, NRHS;
      Complex            B( LDB, * ), D( * ), DL( * ), DU( * );
      // ..

      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      int                J, K;
      Complex            MULT, TEMP, ZDUM;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CGTSV ', -INFO );
         return;
      }

      if (N == 0) return;

      for (K = 1; K <= N - 1; K++) { // 30
         if ( DL( K ) == ZERO ) {

            // Subdiagonal is zero, no elimination is required.

            if ( D( K ) == ZERO ) {

               // Diagonal is zero: set INFO = K and return; a unique
               // solution can not be found.

               INFO = K;
               return;
            }
         } else if ( CABS1( D( K ) ) >= CABS1( DL( K ) ) ) {

            // No row interchange required

            MULT = DL( K ) / D( K );
            D[K+1] = D( K+1 ) - MULT*DU( K );
            for (J = 1; J <= NRHS; J++) { // 10
               B[K+1][J] = B( K+1, J ) - MULT*B( K, J );
            } // 10
            if[K < ( N-1 ) ) DL( K] = ZERO;
         } else {

            // Interchange rows K and K+1

            MULT = D( K ) / DL( K );
            D[K] = DL( K );
            TEMP = D( K+1 );
            D[K+1] = DU( K ) - MULT*TEMP;
            if ( K < ( N-1 ) ) {
               DL[K] = DU( K+1 );
               DU[K+1] = -MULT*DL( K );
            }
            DU[K] = TEMP;
            for (J = 1; J <= NRHS; J++) { // 20
               TEMP = B( K, J );
               B[K][J] = B( K+1, J );
               B[K+1][J] = TEMP - MULT*B( K+1, J );
            } // 20
         }
      } // 30
      if ( D( N ) == ZERO ) {
         INFO = N;
         return;
      }

      // Back solve with the matrix U from the factorization.

      for (J = 1; J <= NRHS; J++) { // 50
         B[N][J] = B( N, J ) / D( N );
         if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
         for (K = N - 2; K >= 1; K--) { // 40
            B[K][J] = ( B( K, J )-DU( K )*B( K+1, J )-DL( K )* B( K+2, J ) ) / D( K );
         } // 40
      } // 50

      }
