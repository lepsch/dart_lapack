      SUBROUTINE ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, TRAN;
      int                I, IASCL, IBSCL, J, MAXMN, BROW, SCLLEN, TSZO, TSZM, LWO, LWM, LW1, LW2, WSIZEO, WSIZEM, INFO2;
      double             ANRM, BIGNUM, BNRM, SMLNUM, DUM( 1 );
      COMPLEX*16         TQ( 5 ), WORKQ( 1 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      // EXTERNAL ZGEQR, ZGEMQR, ZLASCL, ZLASET, ZTRTRS, XERBLA, ZGELQ, ZGEMLQ
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MAXMN = MAX( M, N )
      TRAN  = LSAME( TRANS, 'C' )
*
      LQUERY = ( LWORK.EQ.-1 .OR. LWORK.EQ.-2 )
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
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*     Determine the optimum and minimum LWORK
*
       IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         WSIZEO = 1
         WSIZEM = 1
       ELSE IF( M.GE.N ) THEN
         CALL ZGEQR( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         CALL ZGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 )
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         CALL ZGEQR( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         CALL ZGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 )
         LWM = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       ELSE
         CALL ZGELQ( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         CALL ZGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 )
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         CALL ZGELQ( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         CALL ZGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 )
         LWM  = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       END IF
*
       IF( ( LWORK.LT.WSIZEM ).AND.( .NOT.LQUERY ) ) THEN
          INFO = -10
       END IF
*
       WORK( 1 ) = DBLE( WSIZEO )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGETSLS', -INFO )
        RETURN
      END IF
      IF( LQUERY ) THEN
        IF( LWORK.EQ.-2 ) WORK( 1 ) = DBLE( WSIZEM )
        RETURN
      END IF
      IF( LWORK.LT.WSIZEO ) THEN
        LW1 = TSZM
        LW2 = LWM
      ELSE
        LW1 = TSZO
        LW2 = LWO
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
           CALL ZLASET( 'FULL', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
           RETURN
      END IF
*
*     Get machine parameters
*
       SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
       BIGNUM = ONE / SMLNUM
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', M, N, A, LDA, DUM )
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
         CALL ZLASET( 'F', MAXMN, NRHS, CZERO, CZERO, B, LDB )
         GO TO 50
      END IF
*
      BROW = M
      IF ( TRAN ) THEN
        BROW = N
      END IF
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, DUM )
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
      IF ( M.GE.N ) THEN
*
*        compute QR factorization of A
*
        CALL ZGEQR( M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO )
        IF ( .NOT.TRAN ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
*
          CALL ZGEMQR( 'L' , 'C', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO )
*
*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
          CALL ZTRTRS( 'U', 'N', 'N', N, NRHS, A, LDA, B, LDB, INFO )
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
          SCLLEN = N
        ELSE
*
*           Overdetermined system of equations A**T * X = B
*
*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
*
            CALL ZTRTRS( 'U', 'C', 'N', N, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           B(N+1:M,1:NRHS) = CZERO
*
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = CZERO
   10          CONTINUE
   20       CONTINUE
*
*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
*
            CALL ZGEMQR( 'L', 'N', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO )
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        Compute LQ factorization of A
*
         CALL ZGELQ( M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TRAN ) THEN
*
*           underdetermined system of equations A * X = B
*
*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            CALL ZTRTRS( 'L', 'N', 'N', M, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           B(M+1:N,1:NRHS) = 0
*
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = CZERO
   30          CONTINUE
   40       CONTINUE
*
*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
*
            CALL ZGEMLQ( 'L', 'C', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A**T * X - B ||
*
*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
*
            CALL ZGEMLQ( 'L', 'N', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
*
            CALL ZTRTRS( 'L', 'C', 'N', M, NRHS, A, LDA, B, LDB, INFO )
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
   50 CONTINUE
      WORK( 1 ) = DBLE( TSZO + LWO )
      RETURN
*
*     End of ZGETSLS
*
      END
