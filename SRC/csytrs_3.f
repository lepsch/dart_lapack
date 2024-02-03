      SUBROUTINE CSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )
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
      COMPLEX            A( LDA, * ), B( LDB, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0,0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                I, J, K, KP
      COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSCAL, CSWAP, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
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
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYTRS_3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Begin Upper
*
*        Solve A*X = B, where A = U*D*U**T.
*
*        P**T * B
*
*        Interchange rows K and IPIV(K) of matrix B in the same order
*        that the formation order of IPIV(I) vector for Upper case.
*
*        (We can do the simple loop over IPIV with decrement -1,
*        since the ABS value of IPIV(I) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         DO K = N, 1, -1
            KP = ABS( IPIV( K ) )
            IF( KP.NE.K ) THEN
               CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
         END DO
*
*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
*
         CALL CTRSM( 'L', 'U', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB )
*
*        Compute D \ B -> B   [ D \ (U \P**T * B) ]
*
         I = N
         DO WHILE ( I.GE.1 )
            IF( IPIV( I ).GT.0 ) THEN
               CALL CSCAL( NRHS, ONE / A( I, I ), B( I, 1 ), LDB )
            ELSE IF ( I.GT.1 ) THEN
               AKM1K = E( I )
               AKM1 = A( I-1, I-1 ) / AKM1K
               AK = A( I, I ) / AKM1K
               DENOM = AKM1*AK - ONE
               DO J = 1, NRHS
                  BKM1 = B( I-1, J ) / AKM1K
                  BK = B( I, J ) / AKM1K
                  B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
               END DO
               I = I - 1
            END IF
            I = I - 1
         END DO
*
*        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
*
         CALL CTRSM( 'L', 'U', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB )
*
*        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
*
*        Interchange rows K and IPIV(K) of matrix B in reverse order
*        from the formation order of IPIV(I) vector for Upper case.
*
*        (We can do the simple loop over IPIV with increment 1,
*        since the ABS value of IPIV( I ) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         DO K = 1, N, 1
            KP = ABS( IPIV( K ) )
            IF( KP.NE.K ) THEN
               CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
         END DO
*
      ELSE
*
*        Begin Lower
*
*        Solve A*X = B, where A = L*D*L**T.
*
*        P**T * B
*        Interchange rows K and IPIV(K) of matrix B in the same order
*        that the formation order of IPIV(I) vector for Lower case.
*
*        (We can do the simple loop over IPIV with increment 1,
*        since the ABS value of IPIV(I) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         DO K = 1, N, 1
            KP = ABS( IPIV( K ) )
            IF( KP.NE.K ) THEN
               CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
         END DO
*
*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
*
         CALL CTRSM( 'L', 'L', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB )
*
*        Compute D \ B -> B   [ D \ (L \P**T * B) ]
*
         I = 1
         DO WHILE ( I.LE.N )
            IF( IPIV( I ).GT.0 ) THEN
               CALL CSCAL( NRHS, ONE / A( I, I ), B( I, 1 ), LDB )
            ELSE IF( I.LT.N ) THEN
               AKM1K = E( I )
               AKM1 = A( I, I ) / AKM1K
               AK = A( I+1, I+1 ) / AKM1K
               DENOM = AKM1*AK - ONE
               DO  J = 1, NRHS
                  BKM1 = B( I, J ) / AKM1K
                  BK = B( I+1, J ) / AKM1K
                  B( I, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
               END DO
               I = I + 1
            END IF
            I = I + 1
         END DO
*
*        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
*
         CALL CTRSM('L', 'L', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB )
*
*        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
*
*        Interchange rows K and IPIV(K) of matrix B in reverse order
*        from the formation order of IPIV(I) vector for Lower case.
*
*        (We can do the simple loop over IPIV with decrement -1,
*        since the ABS value of IPIV(I) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         DO K = N, 1, -1
            KP = ABS( IPIV( K ) )
            IF( KP.NE.K ) THEN
               CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
         END DO
*
*        END Lower
*
      END IF
*
      RETURN
*
*     End of CSYTRS_3
*
      END
