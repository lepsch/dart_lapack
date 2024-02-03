      SUBROUTINE SSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, LDA, N
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      REAL               A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                K, KP, KSTEP
      REAL               AK, AKKP1, AKP1, D, T, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SSWAP, SSYMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N ) GO TO 40
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL SCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL SSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( K-1, WORK, 1, A( 1, K ), 1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL SCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL SSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( K-1, WORK, 1, A( 1, K ), 1 )                A( K, K+1 ) = A( K, K+1 ) - SDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL SCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL SSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )                A( K+1, K+1 ) = A( K+1, K+1 ) - SDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            CALL SSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL SSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   50    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 ) GO TO 60
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL SCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL SSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( N-K, WORK, 1, A( K+1, K ), 1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL SCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL SSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( N-K, WORK, 1, A( K+1, K ), 1 )                A( K, K-1 ) = A( K, K-1 ) - SDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 )
               CALL SCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL SSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1 ) = A( K-1, K-1 ) - SDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            IF( KP.LT.N ) CALL SSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL SSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of SSYTRI
*
      END
