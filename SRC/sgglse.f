      SUBROUTINE SGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( * ), D( * ), WORK( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;
*     ..
*     .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGGRQF, SORMQR, SORMRQ, STRMV, STRTRS, XERBLA
*     ..
*     .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.GT.N .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -7
      END IF
*
*     Calculate workspace
*
      IF( INFO.EQ.0) THEN
         IF( N.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            NB1 = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'SGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'SORMQR', ' ', M, N, P, -1 )
            NB4 = ILAENV( 1, 'SORMRQ', ' ', M, N, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = P + MN + MAX( M, N )*NB
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGLSE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Compute the GRQ factorization of matrices B and A:
*
*            B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
*                        N-P  P                     (  0  R22 ) M+P-N
*                                                      N-P  P
*
*     where T12 and R11 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL SGGRQF( P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ), WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = INT( WORK( P+MN+1 ) )
*
*     Update c = Z**T *c = ( c1 ) N-P
*                          ( c2 ) M+P-N
*
      CALL SORMQR( 'Left', 'Transpose', M, 1, MN, A, LDA, WORK( P+1 ), C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
*     Solve T12*x2 = d for x2
*
      IF( P.GT.0 ) THEN
         CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', P, 1, B( 1, N-P+1 ), LDB, D, P, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Put the solution in X
*
         CALL SCOPY( P, D, 1, X( N-P+1 ), 1 )
*
*        Update c1
*
         CALL SGEMV( 'No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D, 1, ONE, C, 1 )
      END IF
*
*     Solve R11*x1 = c1 for x1
*
      IF( N.GT.P ) THEN
         CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA, C, N-P, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
*        Put the solutions in X
*
         CALL SCOPY( N-P, C, 1, X, 1 )
      END IF
*
*     Compute the residual vector:
*
      IF( M.LT.N ) THEN
         NR = M + P - N
         IF( NR.GT.0 ) CALL SGEMV( 'No transpose', NR, N-M, -ONE, A( N-P+1, M+1 ), LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 )
      ELSE
         NR = P
      END IF
      IF( NR.GT.0 ) THEN
         CALL STRMV( 'Upper', 'No transpose', 'Non unit', NR, A( N-P+1, N-P+1 ), LDA, D, 1 )
         CALL SAXPY( NR, -ONE, D, 1, C( N-P+1 ), 1 )
      END IF
*
*     Backward transformation x = Q**T*x
*
      CALL SORMRQ( 'Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ), X, N, WORK( P+MN+1 ), LWORK-P-MN, INFO )
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
      RETURN
*
*     End of SGGLSE
*
      END
