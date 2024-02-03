      SUBROUTINE STFTTP( TRANSR, UPLO, N, ARF, AP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( 0: * ), ARF( 0: * )
      // ..
*
*  =====================================================================
*
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
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STFTTP', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         IF( NORMALTRANSR ) THEN
            AP( 0 ) = ARF( 0 )
         ELSE
            AP( 0 ) = ARF( 0 )
         END IF
         RETURN
      END IF
*
      // Size of array ARF(0:NT-1)
*
      NT = N*( N+1 ) / 2
*
      // Set N1 and N2 depending on LOWER
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
      // If N is odd, set NISODD = .TRUE.
      // If N is even, set K = N/2 and NISODD = .FALSE.
*
      // set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe)
      // where noe = 0 if n is even, noe = 1 if n is odd
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
         LDA = N + 1
      ELSE
         NISODD = .TRUE.
         LDA = N
      END IF
*
      // ARF^C has lda rows and n+1-noe cols
*
      IF( .NOT.NORMALTRANSR ) LDA = ( N+1 ) / 2
*
      // start execution: there are eight cases
*
      IF( NISODD ) THEN
*
         // N is odd
*
         IF( NORMALTRANSR ) THEN
*
            // N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
              // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
              // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
              // T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n
*
               IJP = 0
               JP = 0
               DO J = 0, N2
                  DO I = J, N - 1
                     IJ = I + JP
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JP = JP + LDA
               END DO
               DO I = 0, N2 - 1
                  DO J = 1 + I, N2
                     IJ = I + J*LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
*
            ELSE
*
              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)
*
               IJP = 0
               DO J = 0, N1 - 1
                  IJ = N2 + J
                  DO I = 0, J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                     IJ = IJ + LDA
                  END DO
               END DO
               JS = 0
               DO J = N1, N - 1
                  IJ = JS
                  DO IJ = JS, JS + J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
*
            END IF
*
         ELSE
*
            // N is odd and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
               // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
*
               IJP = 0
               DO I = 0, N2
                  DO IJ = I*( LDA+1 ), N*LDA - 1, LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
               JS = 1
               DO J = 0, N2 - 1
                  DO IJ = JS, JS + N2 - J - 1
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA + 1
               END DO
*
            ELSE
*
               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
               // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
*
               IJP = 0
               JS = N2*LDA
               DO J = 0, N1 - 1
                  DO IJ = JS, JS + J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
               DO I = 0, N1
                  DO IJ = I, I + ( N1+I )*LDA, LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
*
            END IF
*
         END IF
*
      ELSE
*
         // N is even
*
         IF( NORMALTRANSR ) THEN
*
            // N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
               // T1 -> a(1), T2 -> a(0), S -> a(k+1)
*
               IJP = 0
               JP = 0
               DO J = 0, K - 1
                  DO I = J, N - 1
                     IJ = 1 + I + JP
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JP = JP + LDA
               END DO
               DO I = 0, K - 1
                  DO J = I, K - 1
                     IJ = I + J*LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
*
            ELSE
*
               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)
*
               IJP = 0
               DO J = 0, K - 1
                  IJ = K + 1 + J
                  DO I = 0, J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                     IJ = IJ + LDA
                  END DO
               END DO
               JS = 0
               DO J = K, N - 1
                  IJ = JS
                  DO IJ = JS, JS + J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
*
            END IF
*
         ELSE
*
            // N is even and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
*
               IJP = 0
               DO I = 0, K - 1
                  DO IJ = I + ( I+1 )*LDA, ( N+1 )*LDA - 1, LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
               JS = 0
               DO J = 0, K - 1
                  DO IJ = JS, JS + K - J - 1
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA + 1
               END DO
*
            ELSE
*
               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
*
               IJP = 0
               JS = ( K+1 )*LDA
               DO J = 0, K - 1
                  DO IJ = JS, JS + J
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
                  JS = JS + LDA
               END DO
               DO I = 0, K - 1
                  DO IJ = I, I + ( K+I )*LDA, LDA
                     AP( IJP ) = ARF( IJ )
                     IJP = IJP + 1
                  END DO
               END DO
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
      // End of STFTTP
*
      END
