      SUBROUTINE CTFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N, LDA;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( 0: LDA-1, 0: * ), ARF( 0: * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K, NT, NX2, NP1X2;
      int                I, J, L, IJ;
*     ..
*     .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTFTTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) THEN
            IF( NORMALTRANSR ) THEN
               A( 0, 0 ) = ARF( 0 )
            ELSE
               A( 0, 0 ) = CONJG( ARF( 0 ) )
            END IF
         END IF
         RETURN
      END IF
*
*     Size of array ARF(1:2,0:nt-1)
*
      NT = N*( N+1 ) / 2
*
*     set N1 and N2 depending on LOWER: for N even N1=N2=K
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
*     N--by--(N+1)/2.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
         IF( .NOT.LOWER ) NP1X2 = N + N + 2
      ELSE
         NISODD = .TRUE.
         IF( .NOT.LOWER ) NX2 = N + N
      END IF
*
      IF( NISODD ) THEN
*
*        N is odd
*
         IF( NORMALTRANSR ) THEN
*
*           N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n
*
               IJ = 0
               DO J = 0, N2
                  DO I = N1, N2 + J
                     A( N2+J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
*             T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n
*
               IJ = NT - N
               DO J = N - 1, N1, -1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = J - N1, N1 - 1
                     A( J-N1, L ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NX2
               END DO
*
            END IF
*
         ELSE
*
*           N is odd and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is odd
*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
*              T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1
*
               IJ = 0
               DO J = 0, N2 - 1
                  DO I = 0, J
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  DO I = N1 + J, N - 1
                     A( I, N1+J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = N2, N - 1
                  DO I = 0, N1 - 1
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is odd
*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
*              T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda = n2
*
               IJ = 0
               DO J = 0, N1
                  DO I = N1, N - 1
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, N1 - 1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = N2 + J, N - 1
                     A( N2+J, L ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            END IF
*
         END IF
*
      ELSE
*
*        N is even
*
         IF( NORMALTRANSR ) THEN
*
*           N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
*              T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1
*
               IJ = 0
               DO J = 0, K - 1
                  DO I = K, K + J
                     A( K+J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
*              T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1
*
               IJ = NT - N - 1
               DO J = N - 1, K, -1
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = J - K, K - 1
                     A( J-K, L ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NP1X2
               END DO
*
            END IF
*
         ELSE
*
*           N is even and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B)
*              T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) :
*              T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k
*
               IJ = 0
               J = K
               DO I = K, N - 1
                  A( I, J ) = ARF( IJ )
                  IJ = IJ + 1
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
                  DO I = K + 1 + J, N - 1
                     A( I, K+1+J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = K - 1, N - 1
                  DO I = 0, K - 1
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B)
*              T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0)
*              T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k
*
               IJ = 0
               DO J = 0, K
                  DO I = K, N - 1
                     A( J, I ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     A( I, J ) = ARF( IJ )
                     IJ = IJ + 1
                  END DO
                  DO L = K + 1 + J, N - 1
                     A( K+1+J, L ) = CONJG( ARF( IJ ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
*              Note that here J = K-1
*
               DO I = 0, J
                  A( I, J ) = ARF( IJ )
                  IJ = IJ + 1
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
*     End of CTFTTR
*
      END
