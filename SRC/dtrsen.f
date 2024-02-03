      SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI, M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LIWORK, LWORK, M, N;
      double             S, SEP;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double             Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), WR( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS, WANTSP;
      int                IERR, K, KASE, KK, KS, LIWMIN, LWMIN, N1, N2, NN;
      double             EST, RNORM, SCALE;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLANGE;
      // EXTERNAL LSAME, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLACPY, DTREXC, DTRSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
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
      } else {

         // Set M to the dimension of the specified invariant subspace,
         // and test LWORK and LIWORK.

         M = 0
         PAIR = .FALSE.
         DO 10 K = 1, N
            IF( PAIR ) THEN
               PAIR = .FALSE.
            } else {
               IF( K.LT.N ) THEN
                  IF( T( K+1, K ).EQ.ZERO ) THEN
                     IF( SELECT( K ) ) M = M + 1
                  } else {
                     PAIR = .TRUE.
                     IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
                  END IF
               } else {
                  IF( SELECT( N ) ) M = M + 1
               END IF
            END IF
   10    CONTINUE

         N1 = M
         N2 = N - M
         NN = N1*N2

         IF( WANTSP ) THEN
            LWMIN = MAX( 1, 2*NN )
            LIWMIN = MAX( 1, NN )
         ELSE IF( LSAME( JOB, 'N' ) ) THEN
            LWMIN = MAX( 1, N )
            LIWMIN = 1
         ELSE IF( LSAME( JOB, 'E' ) ) THEN
            LWMIN = MAX( 1, NN )
            LIWMIN = 1
         END IF

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -17
         END IF
      END IF

      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible.

      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS ) S = ONE          IF( WANTSP ) SEP = DLANGE( '1', N, N, T, LDT, WORK )
         GO TO 40
      END IF

      // Collect the selected blocks at the top-left corner of T.

      KS = 0
      PAIR = .FALSE.
      DO 20 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         } else {
            SWAP = SELECT( K )
            IF( K.LT.N ) THEN
               IF( T( K+1, K ).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT( K+1 )
               END IF
            END IF
            IF( SWAP ) THEN
               KS = KS + 1

               // Swap the K-th block to position KS.

               IERR = 0
               KK = K
               IF( K.NE.KS ) CALL DTREXC( COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK, IERR )
               IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN

                  // Blocks too close to swap: exit.

                  INFO = 1
                  IF( WANTS ) S = ZERO                   IF( WANTSP ) SEP = ZERO
                  GO TO 40
               END IF
               IF( PAIR ) KS = KS + 1
            END IF
         END IF
   20 CONTINUE

      IF( WANTS ) THEN

         // Solve Sylvester equation for R:

            // T11*R - R*T22 = scale*T12

         CALL DLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )

         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.

         RNORM = DLANGE( 'F', N1, N2, WORK, N1, WORK )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         } else {
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* SQRT( RNORM ) )
         END IF
      END IF

      IF( WANTSP ) THEN

         // Estimate sep(T11,T22).

         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL DLACN2( NN, WORK( NN+1 ), WORK, IWORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN

               // Solve  T11*R - R*T22 = scale*X.

               CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )
            } else {

               // Solve T11**T*R - R*T22**T = scale*X.

               CALL DTRSYL( 'T', 'T', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR )
            END IF
            GO TO 30
         END IF

         SEP = SCALE / EST
      END IF

   40 CONTINUE

      // Store the output eigenvalues in WR and WI.

      DO 50 K = 1, N
         WR( K ) = T( K, K )
         WI( K ) = ZERO
   50 CONTINUE
      DO 60 K = 1, N - 1
         IF( T( K+1, K ).NE.ZERO ) THEN
            WI( K ) = SQRT( ABS( T( K, K+1 ) ) )* SQRT( ABS( T( K+1, K ) ) )
            WI( K+1 ) = -WI( K )
         END IF
   60 CONTINUE

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DTRSEN

      }
