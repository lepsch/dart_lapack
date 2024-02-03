      SUBROUTINE CGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * )
      // ..
*
*  ===================================================================
*
      // .. Parameters ..
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CGGQRF, CTRTRS, CUNMQR, CUNMRQ, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters
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
      // Calculate workspace
*
      IF( INFO.EQ.0) THEN
         IF( N.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            NB1 = ILAENV( 1, 'CGEQRF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'CGERQF', ' ', N, M, -1, -1 )
            NB3 = ILAENV( 1, 'CUNMQR', ' ', N, M, P, -1 )
            NB4 = ILAENV( 1, 'CUNMRQ', ' ', N, M, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = M + NP + MAX( N, P )*NB
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGGLM', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) THEN
         DO I = 1, M
            X(I) = CZERO
         END DO
         DO I = 1, P
            Y(I) = CZERO
         END DO
         RETURN
      END IF
*
      // Compute the GQR factorization of matrices A and B:
*
           // Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M
                    // (  0  ) N-M                 (  0    T22 ) N-M
                       // M                         M+P-N  N-M
*
      // where R11 and T22 are upper triangular, and Q and Z are
      // unitary.
*
      CALL CGGQRF( N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = INT( WORK( M+NP+1 ) )
*
      // Update left-hand-side vector d = Q**H*d = ( d1 ) M
                                                // ( d2 ) N-M
*
      CALL CUNMQR( 'Left', 'Conjugate transpose', N, 1, M, A, LDA, WORK, D, MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
      // Solve T22*y2 = d2 for y2
*
      IF( N.GT.M ) THEN
         CALL CTRTRS( 'Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
*
         CALL CCOPY( N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 )
      END IF
*
      // Set y1 = 0
*
      DO 10 I = 1, M + P - N
         Y( I ) = CZERO
   10 CONTINUE
*
      // Update d1 = d1 - T12*y2
*
      CALL CGEMV( 'No transpose', M, N-M, -CONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, CONE, D, 1 )
*
      // Solve triangular system: R11*x = d1
*
      IF( M.GT.0 ) THEN
         CALL CTRTRS( 'Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
         // Copy D to X
*
         CALL CCOPY( M, D, 1, X, 1 )
      END IF
*
      // Backward transformation y = Z**H *y
*
      CALL CUNMRQ( 'Left', 'Conjugate transpose', P, 1, NP, B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      WORK( 1 ) = M + NP + MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
      RETURN
*
      // End of CGGGLM
*
      END
