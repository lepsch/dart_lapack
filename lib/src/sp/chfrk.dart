// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void chfrk(final int TRANSR, final int UPLO, final int TRANS, final int N, final int K, final int ALPHA, final Matrix<double> A_, final int LDA, final int BETA, final int C,) {
  final A = A_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               ALPHA, BETA;
      int                K, LDA, N;
      String             TRANS, TRANSR, UPLO;
      Complex            A( LDA, * ), C( * );
      // ..

// =====================================================================

      // ..
      // .. Parameters ..
      double               ONE, ZERO;
      Complex            CZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      const              CZERO = ( 0.0, 0.0 ) ;
      bool               LOWER, NORMALTRANSR, NISODD, NOTRANS;
      int                INFO, NROWA, J, NK, N1, N2;
      Complex            CALPHA, CBETA;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, CMPLX


      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      NOTRANS = lsame( TRANS, 'N' );

      if ( NOTRANS ) {
         NROWA = N;
      } else {
         NROWA = K;
      }

      if ( !NORMALTRANSR && !lsame( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( !NOTRANS && !lsame( TRANS, 'C' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, NROWA ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('CHFRK ', -INFO );
         return;
      }

      // Quick return if possible.

      // The quick return case: ((ALPHA == 0) && (BETA != ZERO)) is not
      // done (it is in CHERK for example) and left in the general case.

      if( ( N == 0 ) || ( ( ( ALPHA == ZERO ) || ( K == 0 ) ) && ( BETA == ONE ) ) )return;

      if ( ( ALPHA == ZERO ) && ( BETA == ZERO ) ) {
         for (J = 1; J <= ( ( N*( N+1 ) ) / 2 ); J++) {
            C[J] = CZERO;
         }
         return;
      }

      CALPHA = CMPLX( ALPHA, ZERO );
      CBETA = CMPLX( BETA, ZERO );

      // C is N-by-N.
      // If N is odd, set NISODD = true , and N1 and N2.
      // If N is even, NISODD = false , and NK.

      if ( (N % 2) == 0 ) {
         NISODD = false;
         NK = N / 2;
      } else {
         NISODD = true;
         if ( LOWER ) {
            N2 = N / 2;
            N1 = N - N2;
         } else {
            N1 = N / 2;
            N2 = N - N1;
         }
      }

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'

                  cherk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  cherk('U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N );
                  cgemm('N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'

                  cherk('L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  cherk('U', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N );
                  cgemm('C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N );

               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  cherk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  cherk('U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N );
                  cgemm('N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, CBETA, C( 1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'

                  cherk('L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  cherk('U', 'C', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N );
                  cgemm('C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, CBETA, C( 1 ), N );

               }

            }

         } else {

            // N is odd, and TRANSR = 'C'

            if ( LOWER ) {

               // N is odd, TRANSR = 'C', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'

                  cherk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  cherk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 );
                  cgemm('N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, CBETA, C( N1*N1+1 ), N1 );

               } else {

                  // N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'

                  cherk('U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  cherk('L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 );
                  cgemm('C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, CBETA, C( N1*N1+1 ), N1 );

               }

            } else {

               // N is odd, TRANSR = 'C', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'

                  cherk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  cherk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  cgemm('N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 );

               } else {

                  // N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'

                  cherk('U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  cherk('L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  cgemm('C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 );

               }

            }

         }

      } else {

         // N is even

         if ( NORMALTRANSR ) {

            // N is even and TRANSR = 'N'

            if ( LOWER ) {

               // N is even, TRANSR = 'N', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'

                  cherk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  cherk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 );
                  cgemm('N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'

                  cherk('L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  cherk('U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 );
                  cgemm('C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 );

               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  cherk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  cherk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 );
                  cgemm('N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( 1 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'

                  cherk('L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  cherk('U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 );
                  cgemm('C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( 1 ), N+1 );

               }

            }

         } else {

            // N is even, and TRANSR = 'C'

            if ( LOWER ) {

               // N is even, TRANSR = 'C', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'

                  cherk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  cherk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK );
                  cgemm('N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK );

               } else {

                  // N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'

                  cherk('U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  cherk('L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK );
                  cgemm('C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK );

               }

            } else {

               // N is even, TRANSR = 'C', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'

                  cherk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  cherk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  cgemm('N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK );

               } else {

                  // N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'

                  cherk('U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  cherk('L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  cgemm('C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK );

               }

            }

         }

      }

      }
