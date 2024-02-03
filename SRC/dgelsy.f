      SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4       double             ANRM, BIGNUM, BNRM, C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, WSIZE;;
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL ILAENV, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEQP3, DLAIC1, DLASCL, DLASET, DORMQR, DORMRZ, DTRSM, DTZRZF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1

      // Test the input arguments.

      INFO = 0
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
      END IF

      // Figure out optimal block size

      IF( INFO.EQ.0 ) THEN
         IF( MN.EQ.0 .OR. NRHS.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'DORMQR', ' ', M, N, NRHS, -1 )
            NB4 = ILAENV( 1, 'DORMRQ', ' ', M, N, NRHS, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = MN + MAX( 2*MN, N + 1, MN + NRHS )
            LWKOPT = MAX( LWKMIN, MN + 2*N + NB*( N + 1 ), 2*MN + NB*NRHS )
         END IF
         WORK( 1 ) = LWKOPT

         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELSY', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( MN.EQ.0 .OR. NRHS.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max entries outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN

         // Matrix all zero. Return zero solution.

         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RANK = 0
         GO TO 70
      END IF

      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF

      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R

      CALL DGEQP3( M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )
      WSIZE = MN + WORK( MN+1 )

      // workspace: MN+2*N+NB*(N+1).
      // Details of Householder rotations stored in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK( ISMIN ) = ONE
      WORK( ISMAX ) = ONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
         RANK = 0
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 70
      ELSE
         RANK = 1
      END IF

   10 CONTINUE
      IF( RANK.LT.MN ) THEN
         I = RANK + 1
         CALL DLAIC1( IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 )          CALL DLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 )

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

      // workspace: 3*MN.

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      IF( RANK.LT.N ) CALL DTZRZF( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ), LWORK-2*MN, INFO )

      // workspace: 2*MN.
      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

      CALL DORMQR( 'Left', 'Transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO )
      WSIZE = MAX( WSIZE, 2*MN+WORK( 2*MN+1 ) )

      // workspace: 2*MN+NB*NRHS.

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, ONE, A, LDA, B, LDB )

      DO 40 J = 1, NRHS
         DO 30 I = RANK + 1, N
            B( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE

      // B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)

      IF( RANK.LT.N ) THEN
         CALL DORMRZ( 'Left', 'Transpose', N, NRHS, RANK, N-RANK, A, LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO )
      END IF

      // workspace: 2*MN+NRHS.

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      DO 60 J = 1, NRHS
         DO 50 I = 1, N
            WORK( JPVT( I ) ) = B( I, J )
   50    CONTINUE
         CALL DCOPY( N, WORK( 1 ), 1, B( 1, J ), 1 )
   60 CONTINUE

      // workspace: N.

      // Undo scaling

      IF( IASCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF

   70 CONTINUE
      WORK( 1 ) = LWKOPT

      RETURN

      // End of DGELSY

      }
