// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void clatsy(final int UPLO, final int N, final Matrix<double> X_, final int LDX, final int ISEED,) {
  final X = X_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDX, N;
      int                ISEED( * );
      Complex            X( LDX, * );
      // ..

      Complex            EYE;
      const              EYE = ( 0.0, 1.0 ) ;
      int                I, J, N5;
      double               ALPHA, ALPHA3, BETA;
      Complex            A, B, C, R;
      // ..
      // .. External Functions ..
      //- COMPLEX            CLARND;
      // EXTERNAL CLARND
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      // Initialize constants

      ALPHA = ( 1.+sqrt( 17. ) ) / 8.;
      BETA = ALPHA - 1. / 1000.;
      ALPHA3 = ALPHA*ALPHA*ALPHA;

      // UPLO = 'U':  Upper triangular storage

      if ( UPLO == 'U' ) {

         // Fill the upper triangle of the matrix with zeros.

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               X[I][J] = 0.0;
            } // 10
         } // 20
         N5 = N / 5;
         N5 = N - 5*N5 + 1;

         for (I = N; I >= N5; I -= 5) { // 30
            A = ALPHA3*CLARND( 5, ISEED );
            B = CLARND( 5, ISEED ) / ALPHA;
            C = A - 2.*B*EYE;
            R = C / BETA;
            X[I][I] = A;
            X[I-2][I] = B;
            X[I-2][I-1] = R;
            X[I-2][I-2] = C;
            X[I-1][I-1] = CLARND( 2, ISEED );
            X[I-3][I-3] = CLARND( 2, ISEED );
            X[I-4][I-4] = CLARND( 2, ISEED );
            if ( ( X( I-3, I-3 ) ).abs() > ( X( I-4, I-4 ) ).abs() ) {
               X[I-4][I-3] = 2.0*X( I-3, I-3 );
            } else {
               X[I-4][I-3] = 2.0*X( I-4, I-4 );
            }
         } // 30

         // Clean-up for N not a multiple of 5.

         I = N5 - 1;
         if ( I > 2 ) {
            A = ALPHA3*CLARND( 5, ISEED );
            B = CLARND( 5, ISEED ) / ALPHA;
            C = A - 2.*B*EYE;
            R = C / BETA;
            X[I][I] = A;
            X[I-2][I] = B;
            X[I-2][I-1] = R;
            X[I-2][I-2] = C;
            X[I-1][I-1] = CLARND( 2, ISEED );
            I = I - 3;
         }
         if ( I > 1 ) {
            X[I][I] = CLARND( 2, ISEED );
            X[I-1][I-1] = CLARND( 2, ISEED );
            if ( ( X( I, I ) ).abs() > ( X( I-1, I-1 ) ).abs() ) {
               X[I-1][I] = 2.0*X( I, I );
            } else {
               X[I-1][I] = 2.0*X( I-1, I-1 );
            }
            I = I - 2;
         } else if ( I == 1 ) {
            X[I][I] = CLARND( 2, ISEED );
            I = I - 1;
         }

      // UPLO = 'L':  Lower triangular storage

      } else {

         // Fill the lower triangle of the matrix with zeros.

         for (J = 1; J <= N; J++) { // 50
            for (I = J; I <= N; I++) { // 40
               X[I][J] = 0.0;
            } // 40
         } // 50
         N5 = N / 5;
         N5 = N5*5;

         for (I = 1; I <= N5; I += 5) { // 60
            A = ALPHA3*CLARND( 5, ISEED );
            B = CLARND( 5, ISEED ) / ALPHA;
            C = A - 2.*B*EYE;
            R = C / BETA;
            X[I][I] = A;
            X[I+2][I] = B;
            X[I+2][I+1] = R;
            X[I+2][I+2] = C;
            X[I+1][I+1] = CLARND( 2, ISEED );
            X[I+3][I+3] = CLARND( 2, ISEED );
            X[I+4][I+4] = CLARND( 2, ISEED );
            if ( ( X( I+3, I+3 ) ).abs() > ( X( I+4, I+4 ) ).abs() ) {
               X[I+4][I+3] = 2.0*X( I+3, I+3 );
            } else {
               X[I+4][I+3] = 2.0*X( I+4, I+4 );
            }
         } // 60

         // Clean-up for N not a multiple of 5.

         I = N5 + 1;
         if ( I < N-1 ) {
            A = ALPHA3*CLARND( 5, ISEED );
            B = CLARND( 5, ISEED ) / ALPHA;
            C = A - 2.*B*EYE;
            R = C / BETA;
            X[I][I] = A;
            X[I+2][I] = B;
            X[I+2][I+1] = R;
            X[I+2][I+2] = C;
            X[I+1][I+1] = CLARND( 2, ISEED );
            I = I + 3;
         }
         if ( I < N ) {
            X[I][I] = CLARND( 2, ISEED );
            X[I+1][I+1] = CLARND( 2, ISEED );
            if ( ( X( I, I ) ).abs() > ( X( I+1, I+1 ) ).abs() ) {
               X[I+1][I] = 2.0*X( I, I );
            } else {
               X[I+1][I] = 2.0*X( I+1, I+1 );
            }
            I = I + 2;
         } else if ( I == N ) {
            X[I][I] = CLARND( 2, ISEED );
            I = I + 1;
         }
      }

      }
