      SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                K, LDA, N;
      String             TRANS, TRANSR, UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( * );
      // ..

*  =====================================================================

      // ..
      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NORMALTRANSR, NISODD, NOTRANS;
      int                INFO, NROWA, J, NK, N1, N2;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DGEMM, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = LSAME( TRANSR, 'N' );
      LOWER = LSAME( UPLO, 'L' );
      NOTRANS = LSAME( TRANS, 'N' );

      if ( NOTRANS ) {
         NROWA = N;
      } else {
         NROWA = K;
      }

      if ( !NORMALTRANSR && !LSAME( TRANSR, 'T' ) ) {
         INFO = -1;
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( !NOTRANS && !LSAME( TRANS, 'T' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 ) {
         INFO = -5;
      } else if ( LDA < MAX( 1, NROWA ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DSFRK ', -INFO );
         RETURN;
      }

      // Quick return if possible.

      // The quick return case: ((ALPHA == 0) && (BETA != ZERO)) is not
      // done (it is in DSYRK for example) and left in the general case.

      IF( ( N == 0 ) || ( ( ( ALPHA == ZERO ) || ( K == 0 ) ) && ( BETA == ONE ) ) )RETURN;

      if ( ( ALPHA == ZERO ) && ( BETA == ZERO ) ) {
         for (J = 1; J <= ( ( N*( N+1 ) ) / 2 ); J++) {
            C( J ) = ZERO;
         }
         RETURN;
      }

      // C is N-by-N.
      // If N is odd, set NISODD = true , and N1 and N2.
      // If N is even, NISODD = false , and NK.

      if ( MOD( N, 2 ) == 0 ) {
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

                  dsyrk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  dsyrk('U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N );
                  dgemm('N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'

                  dsyrk('L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  dsyrk('U', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N );
                  dgemm('T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N );

               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  dsyrk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  dsyrk('U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N );
                  dgemm('N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, BETA, C( 1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'

                  dsyrk('L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  dsyrk('U', 'T', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N );
                  dgemm('T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, BETA, C( 1 ), N );

               }

            }

         } else {

            // N is odd, and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'

                  dsyrk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  dsyrk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 );
                  dgemm('N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, BETA, C( N1*N1+1 ), N1 );

               } else {

                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'

                  dsyrk('U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  dsyrk('L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 );
                  dgemm('T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, BETA, C( N1*N1+1 ), N1 );

               }

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'

                  dsyrk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  dsyrk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  dgemm('N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 );

               } else {

                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'

                  dsyrk('U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  dsyrk('L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  dgemm('T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 );

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

                  dsyrk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  dsyrk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 );
                  dgemm('N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'

                  dsyrk('L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  dsyrk('U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 );
                  dgemm('T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );

               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  dsyrk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  dsyrk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 );
                  dgemm('N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'

                  dsyrk('L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  dsyrk('U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 );
                  dgemm('T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 );

               }

            }

         } else {

            // N is even, and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'

                  dsyrk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  dsyrk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK );
                  dgemm('N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK );

               } else {

                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'

                  dsyrk('U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  dsyrk('L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK );
                  dgemm('T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK );

               }

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'

                  dsyrk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  dsyrk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  dgemm('N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK );

               } else {

                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'

                  dsyrk('U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  dsyrk('L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  dgemm('T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK );

               }

            }

         }

      }

      RETURN;

      // End of DSFRK

      }
