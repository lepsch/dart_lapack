      SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             RWORK( * ), S( * ), SEP( * );
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( LDWORK, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D0+0 )
      // ..
      // .. Local Scalars ..
      bool               SOMCON, WANTBH, WANTS, WANTSP;
      String             NORMIN;
      int                I, IERR, IX, J, K, KASE, KS;
      double             BIGNUM, EPS, EST, LNRM, RNRM, SCALE, SMLNUM, XNORM;
      COMPLEX*16         CDUM, PROD
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      COMPLEX*16         DUMMY( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, DZNRM2;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, IZAMAX, DLAMCH, DZNRM2, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDRSCL, ZLACN2, ZLACPY, ZLATRS, ZTREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..
*
      // Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
      // Set M to the number of eigenpairs for which condition numbers are
     t // o be computed.
*
      IF( SOMCON ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) ) M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
*
      INFO = 0
      IF( .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRSNA', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) ) RETURN
         END IF
         IF( WANTS ) S( 1 ) = ONE          IF( WANTSP ) SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF
*
      // Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
      KS = 1
      DO 50 K = 1, N
*
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) ) GO TO 50
         END IF
*
         IF( WANTS ) THEN
*
            // Compute the reciprocal condition number of the k-th
            // eigenvalue.
*
            PROD = ZDOTC( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
            RNRM = DZNRM2( N, VR( 1, KS ), 1 )
            LNRM = DZNRM2( N, VL( 1, KS ), 1 )
            S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
*
         END IF
*
         IF( WANTSP ) THEN
*
            // Estimate the reciprocal condition number of the k-th
            // eigenvector.
*
            // Copy the matrix T to the array WORK and swap the k-th
            // diagonal element to the (1,1) position.
*
            CALL ZLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            CALL ZTREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, K, 1, IERR )
*
            // Form  C = T22 - lambda*I in WORK(2:N,2:N).
*
            DO 20 I = 2, N
               WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20       CONTINUE
*
            // Estimate a lower bound for the 1-norm of inv(C**H). The 1st
            // and (N+1)th columns of WORK are used to store work vectors.
*
            SEP( KS ) = ZERO
            EST = ZERO
            KASE = 0
            NORMIN = 'N'
   30       CONTINUE
            CALL ZLACN2( N-1, WORK( 1, N+1 ), WORK, EST, KASE, ISAVE )
*
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
                  // Solve C**H*x = scale*b
*
                  CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR )
               ELSE
*
                  // Solve C*x = scale*b
*
                  CALL ZLATRS( 'Upper', 'No transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR )
               END IF
               NORMIN = 'Y'
               IF( SCALE.NE.ONE ) THEN
*
                  // Multiply by 1/SCALE if doing so will not cause
                  // overflow.
*
                  IX = IZAMAX( N-1, WORK, 1 )
                  XNORM = CABS1( WORK( IX, 1 ) )
                  IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 40
                  CALL ZDRSCL( N, SCALE, WORK, 1 )
               END IF
               GO TO 30
            END IF
*
            SEP( KS ) = ONE / MAX( EST, SMLNUM )
         END IF
*
   40    CONTINUE
         KS = KS + 1
   50 CONTINUE
      RETURN
*
      // End of ZTRSNA
*
      END
