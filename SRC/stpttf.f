      SUBROUTINE STPTTF( TRANSR, UPLO, N, AP, ARF, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( 0: * ), ARF( 0: * )

*  =====================================================================

      // .. Parameters ..
      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K, NT;
      int                I, J, IJ;
      int                IJP, JP, LDA, JS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) {
         INFO = -1
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      }
      if ( INFO.NE.0 ) {
         xerbla('STPTTF', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         if ( NORMALTRANSR ) {
            ARF( 0 ) = AP( 0 )
         } else {
            ARF( 0 ) = AP( 0 )
         }
         RETURN
      }

      // Size of array ARF(0:NT-1)

      NT = N*( N+1 ) / 2

      // Set N1 and N2 depending on LOWER

      if ( LOWER ) {
         N2 = N / 2
         N1 = N - N2
      } else {
         N1 = N / 2
         N2 = N - N1
      }

      // If N is odd, set NISODD = .TRUE.
      // If N is even, set K = N/2 and NISODD = .FALSE.

      // set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe)
      // where noe = 0 if n is even, noe = 1 if n is odd

      if ( MOD( N, 2 ).EQ.0 ) {
         K = N / 2
         NISODD = .FALSE.
         LDA = N + 1
      } else {
         NISODD = .TRUE.
         LDA = N
      }

      // ARF^C has lda rows and n+1-noe cols

      IF( .NOT.NORMALTRANSR ) LDA = ( N+1 ) / 2

      // start execution: there are eight cases

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               IJP = 0
               JP = 0
               DO J = 0, N2
                  DO I = J, N - 1
                     IJ = I + JP
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JP = JP + LDA
               END DO
               DO I = 0, N2 - 1
                  DO J = 1 + I, N2
                     IJ = I + J*LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               IJP = 0
               DO J = 0, N1 - 1
                  IJ = N2 + J
                  DO I = 0, J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                     IJ = IJ + LDA
                  END DO
               END DO
               JS = 0
               DO J = N1, N - 1
                  IJ = JS
                  DO IJ = JS, JS + J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               IJP = 0
               DO I = 0, N2
                  DO IJ = I*( LDA+1 ), N*LDA - 1, LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO
               JS = 1
               DO J = 0, N2 - 1
                  DO IJ = JS, JS + N2 - J - 1
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA + 1
               END DO

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               IJP = 0
               JS = N2*LDA
               DO J = 0, N1 - 1
                  DO IJ = JS, JS + J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
               DO I = 0, N1
                  DO IJ = I, I + ( N1+I )*LDA, LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO

            }

         }

      } else {

         // N is even

         if ( NORMALTRANSR ) {

            // N is even and TRANSR = 'N'

            if ( LOWER ) {

               // N is even, TRANSR = 'N', and UPLO = 'L'

               IJP = 0
               JP = 0
               DO J = 0, K - 1
                  DO I = J, N - 1
                     IJ = 1 + I + JP
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JP = JP + LDA
               END DO
               DO I = 0, K - 1
                  DO J = I, K - 1
                     IJ = I + J*LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               IJP = 0
               DO J = 0, K - 1
                  IJ = K + 1 + J
                  DO I = 0, J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                     IJ = IJ + LDA
                  END DO
               END DO
               JS = 0
               DO J = K, N - 1
                  IJ = JS
                  DO IJ = JS, JS + J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO

            }

         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               IJP = 0
               DO I = 0, K - 1
                  DO IJ = I + ( I+1 )*LDA, ( N+1 )*LDA - 1, LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO
               JS = 0
               DO J = 0, K - 1
                  DO IJ = JS, JS + K - J - 1
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA + 1
               END DO

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               IJP = 0
               JS = ( K+1 )*LDA
               DO J = 0, K - 1
                  DO IJ = JS, JS + J
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
               DO I = 0, K - 1
                  DO IJ = I, I + ( K+I )*LDA, LDA
                     ARF( IJ ) = AP( IJP )
                     IJP = IJP + 1
                  END DO
               END DO

            }

         }

      }

      RETURN

      // End of STPTTF

      }
