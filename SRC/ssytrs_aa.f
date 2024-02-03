      SUBROUTINE SSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LDB, LWORK, INFO
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      int                K, KP, LWKMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGTSV, SSWAP, SLACPY, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( MIN( N, NRHS ).EQ.0 ) THEN
         LWKMIN = 1
      ELSE
         LWKMIN = 3*N-2
      END IF
*
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
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRS_AA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = SROUNDUP_LWORK( LWKMIN )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( N, NRHS ).EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U**T*T*U.
*
*        1) Forward substitution with U**T
*
         IF( N.GT.1 ) THEN
*
*           Pivot, P**T * B -> B
*
            K = 1
            DO WHILE ( K.LE.N )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K + 1
            END DO
*
*           Compute U**T \ B -> B    [ (U**T \P**T * B) ]
*
            CALL STRSM( 'L', 'U', 'T', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB)
         END IF
*
*        2) Solve with triangular matrix T
*
*        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
*
         CALL SLACPY( 'F', 1, N, A(1, 1), LDA+1, WORK(N), 1)
         IF( N.GT.1 ) THEN
             CALL SLACPY( 'F', 1, N-1, A(1, 2), LDA+1, WORK(1), 1)
             CALL SLACPY( 'F', 1, N-1, A(1, 2), LDA+1, WORK(2*N), 1)
         END IF
         CALL SGTSV(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO)
*
*        3) Backward substitution with U
*
         IF( N.GT.1 ) THEN
*
*
*           Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]
*
            CALL STRSM( 'L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B(2, 1), LDB)
*
*           Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
*
            K = N
            DO WHILE ( K.GE.1 )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K - 1
            END DO
         END IF
*
      ELSE
*
*        Solve A*X = B, where A = L*T*L**T.
*
*        1) Forward substitution with L
*
         IF( N.GT.1 ) THEN
*
*           Pivot, P**T * B -> B
*
            K = 1
            DO WHILE ( K.LE.N )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K + 1
            END DO
*
*           Compute L \ B -> B    [ (L \P**T * B) ]
*
            CALL STRSM( 'L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1), LDA, B(2, 1), LDB)
         END IF
*
*        2) Solve with triangular matrix T
*
*        Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL SLACPY( 'F', 1, N, A(1, 1), LDA+1, WORK(N), 1)
         IF( N.GT.1 ) THEN
             CALL SLACPY( 'F', 1, N-1, A(2, 1), LDA+1, WORK(1), 1)
             CALL SLACPY( 'F', 1, N-1, A(2, 1), LDA+1, WORK(2*N), 1)
         END IF
         CALL SGTSV(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO)
*
*        3) Backward substitution with L**T
*
         IF( N.GT.1 ) THEN
*
*           Compute L**T \ B -> B   [ L**T \ (T \ (L \P**T * B) ) ]
*
            CALL STRSM( 'L', 'L', 'T', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB)
*
*           Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
*
            K = N
            DO WHILE ( K.GE.1 )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K - 1
            END DO
         END IF
*
      END IF
*
      RETURN
*
*     End of SSYTRS_AA
*
      END
