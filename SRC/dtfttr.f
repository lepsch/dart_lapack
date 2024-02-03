      SUBROUTINE DTFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      double             A( 0: LDA-1, 0: * ), ARF( 0: * );
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K, NT, NX2, NP1X2;
      int                I, J, L, IJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
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
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('DTFTTR', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.LE.1 ) {
         if ( N.EQ.1 ) {
            A( 0, 0 ) = ARF( 0 )
         }
         RETURN
      }

      // Size of array ARF(0:nt-1)

      NT = N*( N+1 ) / 2

      // set N1 and N2 depending on LOWER: for N even N1=N2=K

      if ( LOWER ) {
         N2 = N / 2
         N1 = N - N2
      } else {
         N1 = N / 2
         N2 = N - N1
      }

      // If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
      // If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
      // N--by--(N+1)/2.

      if ( MOD( N, 2 ).EQ.0 ) {
         K = N / 2
         NISODD = .FALSE.
         IF( .NOT.LOWER ) NP1X2 = N + N + 2
      } else {
         NISODD = .TRUE.
         IF( .NOT.LOWER ) NX2 = N + N
      }

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               IJ = 0
               DO J = 0, N2
                  DO I = N1, N2 + J
                     A( N2+J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO

            } else {

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N
               DO J = N - 1, N1, -1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = J - N1, N1 - 1
                     A( J-N1, L ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NX2
               END DO

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               IJ = 0
               DO J = 0, N2 - 1
                  DO I = 0, J
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO I = N1 + J, N - 1
                     A( I, N1+J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = N2, N - 1
                  DO I = 0, N1 - 1
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO

            } else {

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               IJ = 0
               DO J = 0, N1
                  DO I = N1, N - 1
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, N1 - 1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = N2 + J, N - 1
                     A( N2+J, L ) = ARF( IJ )
                     IJ = IJ + 1
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

               IJ = 0
               DO J = 0, K - 1
                  DO I = K, K + J
                     A( K+J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO

            } else {

               // N is even, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N - 1
               DO J = N - 1, K, -1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = J - K, K - 1
                     A( J-K, L ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NP1X2
               END DO

            }

         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // N is even, TRANSR = 'T', and UPLO = 'L'

               IJ = 0
               J = K
               DO I = K, N - 1
                  A( I, J ) = ARF( IJ )
                  IJ = IJ + 1
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO I = K + 1 + J, N - 1
                     A( I, K+1+J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = K - 1, N - 1
                  DO I = 0, K - 1
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO

            } else {

               // N is even, TRANSR = 'T', and UPLO = 'U'

               IJ = 0
               DO J = 0, K
                  DO I = K, N - 1
                     A( J, I ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = K + 1 + J, N - 1
                     A( K+1+J, L ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               // Note that here, on exit of the loop, J = K-1
               DO I = 0, J
                  A( I, J ) = ARF( IJ )
                  IJ = IJ + 1
               END DO

            }

         }

      }

      RETURN

      // End of DTFTTR

      }
