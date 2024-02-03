      SUBROUTINE STGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, NCYCLE, P
      REAL               TOLA, TOLB
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                MAXIT
      PARAMETER          ( MAXIT = 40 )
      REAL               ZERO, ONE, HUGENUM
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
*
      bool               INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
      int                I, J, KCYCLE
      REAL               A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, ERROR, GAMMA, RWK, SNQ, SNU, SNV, SSMIN
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAGS2, SLAPLL, SLARTG, SLASET, SROT, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, HUGE
      PARAMETER          ( HUGENUM = HUGE(ZERO) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INITU = LSAME( JOBU, 'I' )
      WANTU = INITU .OR. LSAME( JOBU, 'U' )
*
      INITV = LSAME( JOBV, 'I' )
      WANTV = INITV .OR. LSAME( JOBV, 'V' )
*
      INITQ = LSAME( JOBQ, 'I' )
      WANTQ = INITQ .OR. LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( INITU .OR. WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( INITV .OR. WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( INITQ .OR. WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -18
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -20
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STGSJA', -INFO )
         RETURN
      END IF
*
*     Initialize U, V and Q, if necessary
*
      IF( INITU ) CALL SLASET( 'Full', M, M, ZERO, ONE, U, LDU )       IF( INITV ) CALL SLASET( 'Full', P, P, ZERO, ONE, V, LDV )       IF( INITQ ) CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
*
*     Loop until convergence
*
      UPPER = .FALSE.
      DO 40 KCYCLE = 1, MAXIT
*
         UPPER = .NOT.UPPER
*
         DO 20 I = 1, L - 1
            DO 10 J = I + 1, L
*
               A1 = ZERO
               A2 = ZERO
               A3 = ZERO
               IF( K+I.LE.M ) A1 = A( K+I, N-L+I )                IF( K+J.LE.M ) A3 = A( K+J, N-L+J )
*
               B1 = B( I, N-L+I )
               B3 = B( J, N-L+J )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M ) A2 = A( K+I, N-L+J )
                  B2 = B( I, N-L+J )
               ELSE
                  IF( K+J.LE.M ) A2 = A( K+J, N-L+I )
                  B2 = B( J, N-L+I )
               END IF
*
               CALL SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )
*
*              Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A
*
               IF( K+J.LE.M ) CALL SROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ), LDA, CSU, SNU )
*
*              Update I-th and J-th rows of matrix B: V**T *B
*
               CALL SROT( L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB, CSV, SNV )
*
*              Update (N-L+I)-th and (N-L+J)-th columns of matrices
*              A and B: A*Q and B*Q
*
               CALL SROT( MIN( K+L, M ), A( 1, N-L+J ), 1, A( 1, N-L+I ), 1, CSQ, SNQ )
*
               CALL SROT( L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ, SNQ )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M ) A( K+I, N-L+J ) = ZERO
                  B( I, N-L+J ) = ZERO
               ELSE
                  IF( K+J.LE.M ) A( K+J, N-L+I ) = ZERO
                  B( J, N-L+I ) = ZERO
               END IF
*
*              Update orthogonal matrices U, V, Q, if desired.
*
               IF( WANTU .AND. K+J.LE.M ) CALL SROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU, SNU )
*
               IF( WANTV ) CALL SROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV )
*
               IF( WANTQ ) CALL SROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ, SNQ )
*
   10       CONTINUE
   20    CONTINUE
*
         IF( .NOT.UPPER ) THEN
*
*           The matrices A13 and B13 were lower triangular at the start
*           of the cycle, and are now upper triangular.
*
*           Convergence test: test the parallelism of the corresponding
*           rows of A and B.
*
            ERROR = ZERO
            DO 30 I = 1, MIN( L, M-K )
               CALL SCOPY( L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 )
               CALL SCOPY( L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 )
               CALL SLAPLL( L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN )
               ERROR = MAX( ERROR, SSMIN )
   30       CONTINUE
*
            IF( ABS( ERROR ).LE.MIN( TOLA, TOLB ) ) GO TO 50
         END IF
*
*        End of cycle loop
*
   40 CONTINUE
*
*     The algorithm has not converged after MAXIT cycles.
*
      INFO = 1
      GO TO 100
*
   50 CONTINUE
*
*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
*     Compute the generalized singular value pairs (ALPHA, BETA), and
*     set the triangular matrix R to array A.
*
      DO 60 I = 1, K
         ALPHA( I ) = ONE
         BETA( I ) = ZERO
   60 CONTINUE
*
      DO 70 I = 1, MIN( L, M-K )
*
         A1 = A( K+I, N-L+I )
         B1 = B( I, N-L+I )
         GAMMA = B1 / A1
*
         IF( (GAMMA.LE.HUGENUM).AND.(GAMMA.GE.-HUGENUM) ) THEN
*
*           change sign if necessary
*
            IF( GAMMA.LT.ZERO ) THEN
               CALL SSCAL( L-I+1, -ONE, B( I, N-L+I ), LDB )
               IF( WANTV ) CALL SSCAL( P, -ONE, V( 1, I ), 1 )
            END IF
*
            CALL SLARTG( ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ), RWK )
*
            IF( ALPHA( K+I ).GE.BETA( K+I ) ) THEN
               CALL SSCAL( L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ), LDA )
            ELSE
               CALL SSCAL( L-I+1, ONE / BETA( K+I ), B( I, N-L+I ), LDB )                CALL SCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA )
            END IF
*
         ELSE
*
            ALPHA( K+I ) = ZERO
            BETA( K+I ) = ONE
            CALL SCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA )
*
         END IF
*
   70 CONTINUE
*
*     Post-assignment
*
      DO 80 I = M + 1, K + L
         ALPHA( I ) = ZERO
         BETA( I ) = ONE
   80 CONTINUE
*
      IF( K+L.LT.N ) THEN
         DO 90 I = K + L + 1, N
            ALPHA( I ) = ZERO
            BETA( I ) = ZERO
   90    CONTINUE
      END IF
*
  100 CONTINUE
      NCYCLE = KCYCLE
      RETURN
*
*     End of STGSJA
*
      END
