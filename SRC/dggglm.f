      SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * );
*     ..
*
*  ===================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DGGQRF, DORMQR, DORMRQ, DTRTRS, XERBLA
*     ..
*     .. External Functions ..
      int                ILAENV;
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NP = MIN( N, P )
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
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
            NB1 = ILAENV( 1, 'DGEQRF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'DGERQF', ' ', N, M, -1, -1 )
            NB3 = ILAENV( 1, 'DORMQR', ' ', N, M, P, -1 )
            NB4 = ILAENV( 1, 'DORMRQ', ' ', N, M, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = M + NP + MAX( N, P )*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGGLM', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         DO I = 1, M
            X(I) = ZERO
         END DO
         DO I = 1, P
            Y(I) = ZERO
         END DO
         RETURN
      END IF
*
*     Compute the GQR factorization of matrices A and B:
*
*          Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M
*                   (  0  ) N-M                 (  0    T22 ) N-M
*                      M                         M+P-N  N-M
*
*     where R11 and T22 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL DGGQRF( N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = INT( WORK( M+NP+1 ) )
*
*     Update left-hand-side vector d = Q**T*d = ( d1 ) M
*                                               ( d2 ) N-M
*
      CALL DORMQR( 'Left', 'Transpose', N, 1, M, A, LDA, WORK, D, MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
*     Solve T22*y2 = d2 for y2
*
      IF( N.GT.M ) THEN
         CALL DTRTRS( 'Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
*
         CALL DCOPY( N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 )
      END IF
*
*     Set y1 = 0
*
      DO 10 I = 1, M + P - N
         Y( I ) = ZERO
   10 CONTINUE
*
*     Update d1 = d1 - T12*y2
*
      CALL DGEMV( 'No transpose', M, N-M, -ONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, ONE, D, 1 )
*
*     Solve triangular system: R11*x = d1
*
      IF( M.GT.0 ) THEN
         CALL DTRTRS( 'Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
*        Copy D to X
*
         CALL DCOPY( M, D, 1, X, 1 )
      END IF
*
*     Backward transformation y = Z**T *y
*
      CALL DORMRQ( 'Left', 'Transpose', P, 1, NP, B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      WORK( 1 ) = M + NP + MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
      RETURN
*
*     End of DGGGLM
*
      END
