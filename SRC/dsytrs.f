      SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP
      DOUBLE PRECISION   AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U*D*U**T.
*
*        First solve U*D*X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 ) GO TO 30
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in column K of A.
*
            CALL DGER( K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K - 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Interchange rows K-1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K-1 ) CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in columns K-1 and K of A.
*
            CALL DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )             CALL DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ), LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            AKM1K = A( K-1, K )
            AKM1 = A( K-1, K-1 ) / AKM1K
            AK = A( K, K ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            K = K - 2
         END IF
*
         GO TO 10
   30    CONTINUE
*
*        Next solve U**T *X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   40    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N ) GO TO 50
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Multiply by inv(U**T(K)), where U(K) is the transformation
*           stored in column K of A.
*
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB )
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
*           stored in columns K and K+1 of A.
*
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB )             CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
*
*           Interchange rows K and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 2
         END IF
*
         GO TO 40
   50    CONTINUE
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L**T.
*
*        First solve L*D*X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   60    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N ) GO TO 80
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in column K of A.
*
            IF( K.LT.N ) CALL DGER( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K + 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Interchange rows K+1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K+1 ) CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ), LDB, B( K+2, 1 ), LDB )                CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1, B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
*
*           Multiply by the inverse of the diagonal block.
*
            AKM1K = A( K+1, K )
            AKM1 = A( K, K ) / AKM1K
            AK = A( K+1, K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            K = K + 2
         END IF
*
         GO TO 60
   80    CONTINUE
*
*        Next solve L**T *X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   90    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 ) GO TO 100
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Multiply by inv(L**T(K)), where L(K) is the transformation
*           stored in column K of A.
*
            IF( K.LT.N ) CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
*           stored in columns K-1 and K of A.
*
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )                CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ), LDB )
            END IF
*
*           Interchange rows K and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 2
         END IF
*
         GO TO 90
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of DSYTRS
*
      END
