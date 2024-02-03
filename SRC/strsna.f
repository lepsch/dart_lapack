      SUBROUTINE STRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      REAL               S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               PAIR, SOMCON, WANTBH, WANTS, WANTSP;
      int                I, IERR, IFST, ILST, J, K, KASE, KS, N2, NN;
      REAL               BIGNUM, COND, CS, DELTA, DUMM, EPS, EST, LNRM, MU, PROD, PROD1, PROD2, RNRM, SCALE, SMLNUM, SN
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      REAL               DUMMY( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT, SLAMCH, SLAPY2, SNRM2
      // EXTERNAL LSAME, SDOT, SLAMCH, SLAPY2, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLACPY, SLAQTR, STREXC, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH

      SOMCON = LSAME( HOWMNY, 'S' )

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
      ELSE

         // Set M to the number of eigenpairs for which condition numbers
         // are required, and test MM.

         IF( SOMCON ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
               ELSE
                  IF( K.LT.N ) THEN
                     IF( T( K+1, K ).EQ.ZERO ) THEN
                        IF( SELECT( K ) ) M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
                     END IF
                  ELSE
                     IF( SELECT( N ) ) M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF

         IF( MM.LT.M ) THEN
            INFO = -13
         ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRSNA', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) ) RETURN
         END IF
         IF( WANTS ) S( 1 ) = ONE          IF( WANTSP ) SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      KS = 0
      PAIR = .FALSE.
      DO 60 K = 1, N

         // Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.

         IF( PAIR ) THEN
            PAIR = .FALSE.
            GO TO 60
         ELSE
            IF( K.LT.N ) PAIR = T( K+1, K ).NE.ZERO
         END IF

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         IF( SOMCON ) THEN
            IF( PAIR ) THEN
               IF( .NOT.SELECT( K ) .AND. .NOT.SELECT( K+1 ) ) GO TO 60
            ELSE
               IF( .NOT.SELECT( K ) ) GO TO 60
            END IF
         END IF

         KS = KS + 1

         IF( WANTS ) THEN

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            IF( .NOT.PAIR ) THEN

               // Real eigenvalue.

               PROD = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               RNRM = SNRM2( N, VR( 1, KS ), 1 )
               LNRM = SNRM2( N, VL( 1, KS ), 1 )
               S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
            ELSE

               // Complex eigenvalue.

               PROD1 = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               PROD1 = PROD1 + SDOT( N, VR( 1, KS+1 ), 1, VL( 1, KS+1 ), 1 )
               PROD2 = SDOT( N, VL( 1, KS ), 1, VR( 1, KS+1 ), 1 )
               PROD2 = PROD2 - SDOT( N, VL( 1, KS+1 ), 1, VR( 1, KS ), 1 )                RNRM = SLAPY2( SNRM2( N, VR( 1, KS ), 1 ), SNRM2( N, VR( 1, KS+1 ), 1 ) )                LNRM = SLAPY2( SNRM2( N, VL( 1, KS ), 1 ), SNRM2( N, VL( 1, KS+1 ), 1 ) )
               COND = SLAPY2( PROD1, PROD2 ) / ( RNRM*LNRM )
               S( KS ) = COND
               S( KS+1 ) = COND
            END IF
         END IF

         IF( WANTSP ) THEN

            // Estimate the reciprocal condition number of the k-th
            // eigenvector.

            // Copy the matrix T to the array WORK and swap the diagonal
            // block beginning at T(k,k) to the (1,1) position.

            CALL SLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            IFST = K
            ILST = 1
            CALL STREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST, WORK( 1, N+1 ), IERR )

            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN

               // Could not swap because blocks not well separated

               SCALE = ONE
               EST = BIGNUM
            ELSE

               // Reordering successful

               IF( WORK( 2, 1 ).EQ.ZERO ) THEN

                  // Form C = T22 - lambda*I in WORK(2:N,2:N).

                  DO 20 I = 2, N
                     WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20             CONTINUE
                  N2 = 1
                  NN = N - 1
               ELSE

                  // Triangularize the 2 by 2 block by unitary
                 t // ransformation U = [  cs   i*ss ]
                                     // [ i*ss   cs  ].
                  // such that the (1,1) position of WORK is complex
                  // eigenvalue lambda with positive imaginary part. (2,2)
                  // position of WORK is the complex eigenvalue lambda
                  // with negative imaginary  part.

                  MU = SQRT( ABS( WORK( 1, 2 ) ) )* SQRT( ABS( WORK( 2, 1 ) ) )
                  DELTA = SLAPY2( MU, WORK( 2, 1 ) )
                  CS = MU / DELTA
                  SN = -WORK( 2, 1 ) / DELTA

                  // Form

                  // C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
                                           // [   mu                     ]
                                           // [         ..               ]
                                           // [             ..           ]
                                           // [                  mu      ]
                  // where C**T is transpose of matrix C,
                  // and RWORK is stored starting in the N+1-st column of
                  // WORK.

                  DO 30 J = 3, N
                     WORK( 2, J ) = CS*WORK( 2, J )
                     WORK( J, J ) = WORK( J, J ) - WORK( 1, 1 )
   30             CONTINUE
                  WORK( 2, 2 ) = ZERO

                  WORK( 1, N+1 ) = TWO*MU
                  DO 40 I = 2, N - 1
                     WORK( I, N+1 ) = SN*WORK( 1, I+1 )
   40             CONTINUE
                  N2 = 2
                  NN = 2*( N-1 )
               END IF

               // Estimate norm(inv(C**T))

               EST = ZERO
               KASE = 0
   50          CONTINUE
               CALL SLACN2( NN, WORK( 1, N+2 ), WORK( 1, N+4 ), IWORK, EST, KASE, ISAVE )
               IF( KASE.NE.0 ) THEN
                  IF( KASE.EQ.1 ) THEN
                     IF( N2.EQ.1 ) THEN

                        // Real eigenvalue: solve C**T*x = scale*c.

                        CALL SLAQTR( .TRUE., .TRUE., N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     ELSE

                        // Complex eigenvalue: solve
                        // C**T*(p+iq) = scale*(c+id) in real arithmetic.

                        CALL SLAQTR( .TRUE., .FALSE., N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     END IF
                  ELSE
                     IF( N2.EQ.1 ) THEN

                        // Real eigenvalue: solve C*x = scale*c.

                        CALL SLAQTR( .FALSE., .TRUE., N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     ELSE

                        // Complex eigenvalue: solve
                        // C*(p+iq) = scale*(c+id) in real arithmetic.

                        CALL SLAQTR( .FALSE., .FALSE., N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )

                     END IF
                  END IF

                  GO TO 50
               END IF
            END IF

            SEP( KS ) = SCALE / MAX( EST, SMLNUM )
            IF( PAIR ) SEP( KS+1 ) = SEP( KS )
         END IF

         IF( PAIR ) KS = KS + 1

   60 CONTINUE
      RETURN

      // End of STRSNA

      END
