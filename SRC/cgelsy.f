      SUBROUTINE CGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                IMAX, IMIN;
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKOPT, MN, NB, NB1, NB2, NB3, NB4       REAL               ANRM, BIGNUM, BNRM, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, WSIZE;
      COMPLEX            C1, C2, S1, S2
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEQP3, CLAIC1, CLASCL, CLASET, CTRSM, CTZRZF, CUNMQR, CUNMRZ, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, ILAENV, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, CMPLX
      // ..
      // .. Executable Statements ..
*
      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1
*
      // Test the input arguments.
*
      INFO = 0
      NB1 = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )
      NB2 = ILAENV( 1, 'CGERQF', ' ', M, N, -1, -1 )
      NB3 = ILAENV( 1, 'CUNMQR', ' ', M, N, NRHS, -1 )
      NB4 = ILAENV( 1, 'CUNMRQ', ' ', M, N, NRHS, -1 )
      NB = MAX( NB1, NB2, NB3, NB4 )
      LWKOPT = MAX( 1, MN+2*N+NB*(N+1), 2*MN+NB*NRHS )
      WORK( 1 ) = CMPLX( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -7
      ELSE IF( LWORK.LT.( MN+MAX( 2*MN, N+1, MN+NRHS ) ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGELSY', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
*
      // Get machine parameters
*
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
*
      // Scale A, B if max entries outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
         // Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
         // Scale matrix norm down to BIGNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
         // Matrix all zero. Return zero solution.
*
         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         RANK = 0
         GO TO 70
      END IF
*
      BNRM = CLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
         // Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
         // Scale matrix norm down to BIGNUM
*
         CALL CLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R
*
      CALL CGEQP3( M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), LWORK-MN, RWORK, INFO )
      WSIZE = MN + REAL( WORK( MN+1 ) )
*
      // complex workspace: MN+NB*(N+1). real workspace 2*N.
      // Details of Householder rotations stored in WORK(1:MN).
*
      // Determine RANK using incremental condition estimation
*
      WORK( ISMIN ) = CONE
      WORK( ISMAX ) = CONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
         RANK = 0
         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         GO TO 70
      ELSE
         RANK = 1
      END IF
*
   10 CONTINUE
      IF( RANK.LT.MN ) THEN
         I = RANK + 1
         CALL CLAIC1( IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 )          CALL CLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 )
*
         IF( SMAXPR*RCOND.LE.SMINPR ) THEN
            DO 20 I = 1, RANK
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
   20       CONTINUE
            WORK( ISMIN+RANK ) = C1
            WORK( ISMAX+RANK ) = C2
            SMIN = SMINPR
            SMAX = SMAXPR
            RANK = RANK + 1
            GO TO 10
         END IF
      END IF
*
      // complex workspace: 3*MN.
*
      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)
*
      // [R11,R12] = [ T11, 0 ] * Y
*
      IF( RANK.LT.N ) CALL CTZRZF( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ), LWORK-2*MN, INFO )
*
      // complex workspace: 2*MN.
      // Details of Householder rotations stored in WORK(MN+1:2*MN)
*
      // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
*
      CALL CUNMQR( 'Left', 'Conjugate transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO )
      WSIZE = MAX( WSIZE, 2*MN+REAL( WORK( 2*MN+1 ) ) )
*
      // complex workspace: 2*MN+NB*NRHS.
*
      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
*
      CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, CONE, A, LDA, B, LDB )
*
      DO 40 J = 1, NRHS
         DO 30 I = RANK + 1, N
            B( I, J ) = CZERO
   30    CONTINUE
   40 CONTINUE
*
      // B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)
*
      IF( RANK.LT.N ) THEN
         CALL CUNMRZ( 'Left', 'Conjugate transpose', N, NRHS, RANK, N-RANK, A, LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO )
      END IF
*
      // complex workspace: 2*MN+NRHS.
*
      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
*
      DO 60 J = 1, NRHS
         DO 50 I = 1, N
            WORK( JPVT( I ) ) = B( I, J )
   50    CONTINUE
         CALL CCOPY( N, WORK( 1 ), 1, B( 1, J ), 1 )
   60 CONTINUE
*
      // complex workspace: N.
*
      // Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL CLASCL( 'U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL CLASCL( 'U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
*
   70 CONTINUE
      WORK( 1 ) = CMPLX( LWKOPT )
*
      RETURN
*
      // End of CGELSY
*
      END
