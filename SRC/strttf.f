      SUBROUTINE STRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      REAL               A( 0: LDA-1, 0: * ), ARF( 0: * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                I, IJ, J, K, L, N1, N2, NT, NX2, NP1X2;
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
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRTTF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) THEN
            ARF( 0 ) = A( 0, 0 )
         END IF
         RETURN
      END IF

      // Size of array ARF(0:nt-1)

      NT = N*( N+1 ) / 2

      // Set N1 and N2 depending on LOWER: for N even N1=N2=K

      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF

      // If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
      // If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
      // N--by--(N+1)/2.

      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
         IF( .NOT.LOWER ) NP1X2 = N + N + 2
      ELSE
         NISODD = .TRUE.
         IF( .NOT.LOWER ) NX2 = N + N
      END IF

      IF( NISODD ) THEN

         // N is odd

         IF( NORMALTRANSR ) THEN

            // N is odd and TRANSR = 'N'

            IF( LOWER ) THEN

               // N is odd, TRANSR = 'N', and UPLO = 'L'

               IJ = 0
               DO J = 0, N2
                  DO I = N1, N2 + J
                     ARF( IJ ) = A( N2+J, I )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
               END DO

            ELSE

               // N is odd, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N
               DO J = N - 1, N1, -1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = J - N1, N1 - 1
                     ARF( IJ ) = A( J-N1, L )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NX2
               END DO

            END IF

         ELSE

            // N is odd and TRANSR = 'T'

            IF( LOWER ) THEN

               // N is odd, TRANSR = 'T', and UPLO = 'L'

               IJ = 0
               DO J = 0, N2 - 1
                  DO I = 0, J
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
                  DO I = N1 + J, N - 1
                     ARF( IJ ) = A( I, N1+J )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = N2, N - 1
                  DO I = 0, N1 - 1
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
               END DO

            ELSE

               // N is odd, TRANSR = 'T', and UPLO = 'U'

               IJ = 0
               DO J = 0, N1
                  DO I = N1, N - 1
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, N1 - 1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = N2 + J, N - 1
                     ARF( IJ ) = A( N2+J, L )
                     IJ = IJ + 1
                  END DO
               END DO

            END IF

         END IF

      ELSE

         // N is even

         IF( NORMALTRANSR ) THEN

            // N is even and TRANSR = 'N'

            IF( LOWER ) THEN

               // N is even, TRANSR = 'N', and UPLO = 'L'

               IJ = 0
               DO J = 0, K - 1
                  DO I = K, K + J
                     ARF( IJ ) = A( K+J, I )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
               END DO

            ELSE

               // N is even, TRANSR = 'N', and UPLO = 'U'

               IJ = NT - N - 1
               DO J = N - 1, K, -1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = J - K, K - 1
                     ARF( IJ ) = A( J-K, L )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NP1X2
               END DO

            END IF

         ELSE

            // N is even and TRANSR = 'T'

            IF( LOWER ) THEN

               // N is even, TRANSR = 'T', and UPLO = 'L'

               IJ = 0
               J = K
               DO I = K, N - 1
                  ARF( IJ ) = A( I, J )
                  IJ = IJ + 1
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
                  DO I = K + 1 + J, N - 1
                     ARF( IJ ) = A( I, K+1+J )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = K - 1, N - 1
                  DO I = 0, K - 1
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
               END DO

            ELSE

               // N is even, TRANSR = 'T', and UPLO = 'U'

               IJ = 0
               DO J = 0, K
                  DO I = K, N - 1
                     ARF( IJ ) = A( J, I )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = K + 1 + J, N - 1
                     ARF( IJ ) = A( K+1+J, L )
                     IJ = IJ + 1
                  END DO
               END DO
               // Note that here, on exit of the loop, J = K-1
               DO I = 0, J
                  ARF( IJ ) = A( I, J )
                  IJ = IJ + 1
               END DO

            END IF

         END IF

      END IF

      RETURN

      // End of STRTTF

      }
