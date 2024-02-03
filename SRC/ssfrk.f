      SUBROUTINE SSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               ALPHA, BETA
      int                K, LDA, N;
      String             TRANS, TRANSR, UPLO;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), C( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
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
      // EXTERNAL SGEMM, SSYRK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )

      if ( NOTRANS ) {
         NROWA = N
      } else {
         NROWA = K
      }

      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) {
         INFO = -1
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, NROWA ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('SSFRK ', -INFO );
         RETURN
      }

      // Quick return if possible.

      // The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
      // done (it is in SSYRK for example) and left in the general case.

      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND. ( BETA.EQ.ONE ) ) )RETURN

      if ( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ZERO ) ) {
         for (J = 1; J <= ( ( N*( N+1 ) ) / 2 ); J++) {
            C( J ) = ZERO
         }
         RETURN
      }

      // C is N-by-N.
      // If N is odd, set NISODD = true , and N1 and N2.
      // If N is even, NISODD = false , and NK.

      if ( MOD( N, 2 ).EQ.0 ) {
         NISODD = false;
         NK = N / 2
      } else {
         NISODD = true;
         if ( LOWER ) {
            N2 = N / 2
            N1 = N - N2
         } else {
            N1 = N / 2
            N2 = N - N1
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

                  ssyrk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  ssyrk('U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N );
                  sgemm('N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'

                  ssyrk('L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N );
                  ssyrk('U', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N );
                  sgemm('T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N );

               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  ssyrk('L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  ssyrk('U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N );
                  sgemm('N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, BETA, C( 1 ), N );

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'

                  ssyrk('L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N );
                  ssyrk('U', 'T', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N );
                  sgemm('T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, BETA, C( 1 ), N );

               }

            }

         } else {

            // N is odd, and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'

                  ssyrk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  ssyrk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 );
                  sgemm('N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, BETA, C( N1*N1+1 ), N1 );

               } else {

                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'

                  ssyrk('U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 );
                  ssyrk('L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 );
                  sgemm('T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, BETA, C( N1*N1+1 ), N1 );

               }

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'

                  ssyrk('U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  ssyrk('L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  sgemm('N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 );

               } else {

                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'

                  ssyrk('U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 );
                  ssyrk('L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 );
                  sgemm('T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 );

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

                  ssyrk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  ssyrk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 );
                  sgemm('N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'

                  ssyrk('L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 );
                  ssyrk('U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 );
                  sgemm('T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );

               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  ssyrk('L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  ssyrk('U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 );
                  sgemm('N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 );

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'

                  ssyrk('L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 );
                  ssyrk('U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 );
                  sgemm('T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 );

               }

            }

         } else {

            // N is even, and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'

                  ssyrk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  ssyrk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK );
                  sgemm('N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK );

               } else {

                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'

                  ssyrk('U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK );
                  ssyrk('L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK );
                  sgemm('T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK );

               }

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'

                  ssyrk('U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  ssyrk('L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  sgemm('N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK );

               } else {

                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'

                  ssyrk('U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK );
                  ssyrk('L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK );
                  sgemm('T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK );

               }

            }

         }

      }

      RETURN

      // End of SSFRK

      }
