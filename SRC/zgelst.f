      SUBROUTINE ZGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, LWOPT, MN, MNNRHS, NB, NBMIN, SCLLEN
      double             ANRM, BIGNUM, BNRM, SMLNUM;
*     ..
*     .. Local Arrays ..
      double             RWORK( 1 );
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV
      double             DLAMCH, ZLANGE;
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGELQT, ZGEQRT, ZGEMLQT, ZGEMQRT, ZLASCL, ZLASET, ZTRTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
*     Figure out optimal block size and optimal workspace size
*
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
*
         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) ) TPSD = .FALSE.
*
         NB = ILAENV( 1, 'ZGELST', ' ', M, N, -1, -1 )
*
         MNNRHS = MAX( MN, NRHS )
         LWOPT = MAX( 1, (MN+MNNRHS)*NB )
         WORK( 1 ) = DBLE( LWOPT )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGELST ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL ZLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         WORK( 1 ) = DBLE( LWOPT )
         RETURN
      END IF
*
*     *GEQRT and *GELQT routines cannot accept NB larger than min(M,N)
*
      IF( NB.GT.MN ) NB = MN
*
*     Determine the block size from the supplied LWORK
*     ( at this stage we know that LWORK >= (minimum required workspace,
*     but it may be less than optimal)
*
      NB = MIN( NB, LWORK/( MN + MNNRHS ) )
*
*     The minimum value of NB, when blocked code is used
*
      NBMIN = MAX( 2, ILAENV( 2, 'ZGELST', ' ', M, N, -1, -1 ) )
*
      IF( NB.LT.NBMIN ) THEN
         NB = 1
      END IF
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL ZLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         WORK( 1 ) = DBLE( LWOPT )
         RETURN
      END IF
*
      BROW = M
      IF( TPSD ) BROW = N
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL ZLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL ZLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
      IF( M.GE.N ) THEN
*
*        M > N:
*        Compute the blocked QR factorization of A,
*        using the compact WY representation of Q,
*        workspace at least N, optimally N*NB.
*
         CALL ZGEQRT( M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO )
*
         IF( .NOT.TPSD ) THEN
*
*           M > N, A is not transposed:
*           Overdetermined system of equations,
*           least-squares problem, min || A * X - B ||.
*
*           Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS),
*           using the compact WY representation of Q,
*           workspace at least NRHS, optimally NRHS*NB.
*
            CALL ZGEMQRT( 'Left', 'Conjugate transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO )
*
*           Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
            CALL ZTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
            SCLLEN = N
*
         ELSE
*
*           M > N, A is transposed:
*           Underdetermined system of equations,
*           minimum norm solution of A**T * X = B.
*
*           Compute B := inv(R**T) * B in two row blocks of B.
*
*           Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
*
            CALL ZTRTRS( 'Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           Block 2: Zero out all rows below the N-th row in B:
*           B(N+1:M,1:NRHS) = ZERO
*
            DO  J = 1, NRHS
               DO I = N + 1, M
                  B( I, J ) = ZERO
               END DO
            END DO
*
*           Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS),
*           using the compact WY representation of Q,
*           workspace at least NRHS, optimally NRHS*NB.
*
            CALL ZGEMQRT( 'Left', 'No transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO )
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        M < N:
*        Compute the blocked LQ factorization of A,
*        using the compact WY representation of Q,
*        workspace at least M, optimally M*NB.
*
         CALL ZGELQT( M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO )
*
         IF( .NOT.TPSD ) THEN
*
*           M < N, A is not transposed:
*           Underdetermined system of equations,
*           minimum norm solution of A * X = B.
*
*           Compute B := inv(L) * B in two row blocks of B.
*
*           Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            CALL ZTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           Block 2: Zero out all rows below the M-th row in B:
*           B(M+1:N,1:NRHS) = ZERO
*
            DO J = 1, NRHS
               DO I = M + 1, N
                  B( I, J ) = ZERO
               END DO
            END DO
*
*           Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS),
*           using the compact WY representation of Q,
*           workspace at least NRHS, optimally NRHS*NB.
*
            CALL ZGEMLQT( 'Left', 'Conjugate transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO )
*
            SCLLEN = N
*
         ELSE
*
*           M < N, A is transposed:
*           Overdetermined system of equations,
*           least-squares problem, min || A**T * X - B ||.
*
*           Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS),
*           using the compact WY representation of Q,
*           workspace at least NRHS, optimally NRHS*NB.
*
            CALL ZGEMLQT( 'Left', 'No transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1), INFO )
*
*           Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
*
            CALL ZTRTRS( 'Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
            SCLLEN = M
*
         END IF
*
      END IF
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL ZLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL ZLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      END IF
*
      WORK( 1 ) = DBLE( LWOPT )
*
      RETURN
*
*     End of ZGELST
*
      END
