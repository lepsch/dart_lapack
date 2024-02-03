      SUBROUTINE ZHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                K, LDA, N;
      String             TRANS, TRANSR, UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      COMPLEX*16         CZERO
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NORMALTRANSR, NISODD, NOTRANS;
      int                INFO, NROWA, J, NK, N1, N2;
      COMPLEX*16         CALPHA, CBETA
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZHERK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, DCMPLX
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

      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) {
         INFO = -1
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, NROWA ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHFRK ', -INFO )
         RETURN
      }

      // Quick return if possible.

      // The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
      // done (it is in ZHERK for example) and left in the general case.

      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND. ( BETA.EQ.ONE ) ) )RETURN

      if ( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ZERO ) ) {
         DO J = 1, ( ( N*( N+1 ) ) / 2 )
            C( J ) = CZERO
         END DO
         RETURN
      }

      CALPHA = DCMPLX( ALPHA, ZERO )
      CBETA = DCMPLX( BETA, ZERO )

      // C is N-by-N.
      // If N is odd, set NISODD = .TRUE., and N1 and N2.
      // If N is even, NISODD = .FALSE., and NK.

      if ( MOD( N, 2 ).EQ.0 ) {
         NISODD = .FALSE.
         NK = N / 2
      } else {
         NISODD = .TRUE.
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

                  CALL ZHERK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL ZHERK( 'U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N )                   CALL ZGEMM( 'N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N )

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'

                  CALL ZHERK( 'L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL ZHERK( 'U', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N )                   CALL ZGEMM( 'C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N )

               }

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  CALL ZHERK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL ZHERK( 'U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N )                   CALL ZGEMM( 'N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, CBETA, C( 1 ), N )

               } else {

                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'

                  CALL ZHERK( 'L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL ZHERK( 'U', 'C', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N )                   CALL ZGEMM( 'C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, CBETA, C( 1 ), N )

               }

            }

         } else {

            // N is odd, and TRANSR = 'C'

            if ( LOWER ) {

               // N is odd, TRANSR = 'C', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'

                  CALL ZHERK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL ZHERK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 )                   CALL ZGEMM( 'N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, CBETA, C( N1*N1+1 ), N1 )

               } else {

                  // N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'

                  CALL ZHERK( 'U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL ZHERK( 'L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 )                   CALL ZGEMM( 'C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, CBETA, C( N1*N1+1 ), N1 )

               }

            } else {

               // N is odd, TRANSR = 'C', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'

                  CALL ZHERK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL ZHERK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL ZGEMM( 'N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 )

               } else {

                  // N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'

                  CALL ZHERK( 'U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL ZHERK( 'L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL ZGEMM( 'C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 )

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

                  CALL ZHERK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL ZHERK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 )                   CALL ZGEMM( 'N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 )

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'

                  CALL ZHERK( 'L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL ZHERK( 'U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 )                   CALL ZGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 )

               }

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'

                  CALL ZHERK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL ZHERK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL ZGEMM( 'N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( 1 ), N+1 )

               } else {

                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'

                  CALL ZHERK( 'L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL ZHERK( 'U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL ZGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( 1 ), N+1 )

               }

            }

         } else {

            // N is even, and TRANSR = 'C'

            if ( LOWER ) {

               // N is even, TRANSR = 'C', and UPLO = 'L'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'

                  CALL ZHERK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL ZHERK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK )                   CALL ZGEMM( 'N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK )

               } else {

                  // N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'

                  CALL ZHERK( 'U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL ZHERK( 'L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK )                   CALL ZGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK )

               }

            } else {

               // N is even, TRANSR = 'C', and UPLO = 'U'

               if ( NOTRANS ) {

                  // N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'

                  CALL ZHERK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL ZHERK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL ZGEMM( 'N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK )

               } else {

                  // N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'

                  CALL ZHERK( 'U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL ZHERK( 'L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL ZGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK )

               }

            }

         }

      }

      RETURN

      // End of ZHFRK

      }
