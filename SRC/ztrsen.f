      SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S, SEP, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LWORK, M, N;
      double             S, SEP;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      COMPLEX*16         Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTBH, WANTQ, WANTS, WANTSP;
      int                IERR, K, KASE, KS, LWMIN, N1, N2, NN;
      double             EST, RNORM, SCALE;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             ZLANGE;
      // EXTERNAL LSAME, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACN2, ZLACPY, ZTREXC, ZTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..
*
      // Decode and test the input parameters.
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )
*
      // Set M to the number of selected eigenvalues.
*
      M = 0
      DO 10 K = 1, N
         IF( SELECT( K ) ) M = M + 1
   10 CONTINUE
*
      N1 = M
      N2 = N - M
      NN = N1*N2
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( WANTSP ) THEN
         LWMIN = MAX( 1, 2*NN )
      ELSE IF( LSAME( JOB, 'N' ) ) THEN
         LWMIN = 1
      ELSE IF( LSAME( JOB, 'E' ) ) THEN
         LWMIN = MAX( 1, NN )
      END IF
*
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -14
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS ) S = ONE          IF( WANTSP ) SEP = ZLANGE( '1', N, N, T, LDT, RWORK )
         GO TO 40
      END IF
*
      // Collect the selected eigenvalues at the top left corner of T.
*
      KS = 0
      DO 20 K = 1, N
         IF( SELECT( K ) ) THEN
            KS = KS + 1
*
            // Swap the K-th eigenvalue to position KS.
*
            IF( K.NE.KS ) CALL ZTREXC( COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR )
         END IF
   20 CONTINUE
*
      IF( WANTS ) THEN
*
         // Solve the Sylvester equation for R:
*
            // T11*R - R*T22 = scale*T12
*
         CALL ZLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL ZTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )
*
         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.
*
         RNORM = ZLANGE( 'F', N1, N2, WORK, N1, RWORK )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* SQRT( RNORM ) )
         END IF
      END IF
*
      IF( WANTSP ) THEN
*
         // Estimate sep(T11,T22).
*
         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL ZLACN2( NN, WORK( NN+1 ), WORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
               // Solve T11*R - R*T22 = scale*X.
*
               CALL ZTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )
            ELSE
*
               // Solve T11**H*R - R*T22**H = scale*X.
*
               CALL ZTRSYL( 'C', 'C', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )
            END IF
            GO TO 30
         END IF
*
         SEP = SCALE / EST
      END IF
*
   40 CONTINUE
*
      // Copy reordered eigenvalues to W.
*
      DO 50 K = 1, N
         W( K ) = T( K, K )
   50 CONTINUE
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
      // End of ZTRSEN
*
      END
