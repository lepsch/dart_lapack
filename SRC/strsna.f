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
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
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
      if ( .NOT.WANTS .AND. .NOT.WANTSP ) {
         INFO = -1
      } else if ( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) {
         INFO = -8
      } else if ( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) {
         INFO = -10
      } else {

         // Set M to the number of eigenpairs for which condition numbers
         // are required, and test MM.

         if ( SOMCON ) {
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               if ( PAIR ) {
                  PAIR = .FALSE.
               } else {
                  if ( K.LT.N ) {
                     if ( T( K+1, K ).EQ.ZERO ) {
                        IF( SELECT( K ) ) M = M + 1
                     } else {
                        PAIR = .TRUE.
                        IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
                     }
                  } else {
                     IF( SELECT( N ) ) M = M + 1
                  }
               }
   10       CONTINUE
         } else {
            M = N
         }

         if ( MM.LT.M ) {
            INFO = -13
         } else if ( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) {
            INFO = -16
         }
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'STRSNA', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         if ( SOMCON ) {
            IF( .NOT.SELECT( 1 ) ) RETURN
         }
         IF( WANTS ) S( 1 ) = ONE          IF( WANTSP ) SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      KS = 0
      PAIR = .FALSE.
      DO 60 K = 1, N

         // Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.

         if ( PAIR ) {
            PAIR = .FALSE.
            GO TO 60
         } else {
            IF( K.LT.N ) PAIR = T( K+1, K ).NE.ZERO
         }

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         if ( SOMCON ) {
            if ( PAIR ) {
               IF( .NOT.SELECT( K ) .AND. .NOT.SELECT( K+1 ) ) GO TO 60
            } else {
               IF( .NOT.SELECT( K ) ) GO TO 60
            }
         }

         KS = KS + 1

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            if ( .NOT.PAIR ) {

               // Real eigenvalue.

               PROD = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               RNRM = SNRM2( N, VR( 1, KS ), 1 )
               LNRM = SNRM2( N, VL( 1, KS ), 1 )
               S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
            } else {

               // Complex eigenvalue.

               PROD1 = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               PROD1 = PROD1 + SDOT( N, VR( 1, KS+1 ), 1, VL( 1, KS+1 ), 1 )
               PROD2 = SDOT( N, VL( 1, KS ), 1, VR( 1, KS+1 ), 1 )
               PROD2 = PROD2 - SDOT( N, VL( 1, KS+1 ), 1, VR( 1, KS ), 1 )                RNRM = SLAPY2( SNRM2( N, VR( 1, KS ), 1 ), SNRM2( N, VR( 1, KS+1 ), 1 ) )                LNRM = SLAPY2( SNRM2( N, VL( 1, KS ), 1 ), SNRM2( N, VL( 1, KS+1 ), 1 ) )
               COND = SLAPY2( PROD1, PROD2 ) / ( RNRM*LNRM )
               S( KS ) = COND
               S( KS+1 ) = COND
            }
         }

         if ( WANTSP ) {

            // Estimate the reciprocal condition number of the k-th
            // eigenvector.

            // Copy the matrix T to the array WORK and swap the diagonal
            // block beginning at T(k,k) to the (1,1) position.

            CALL SLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            IFST = K
            ILST = 1
            CALL STREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST, WORK( 1, N+1 ), IERR )

            if ( IERR.EQ.1 .OR. IERR.EQ.2 ) {

               // Could not swap because blocks not well separated

               SCALE = ONE
               EST = BIGNUM
            } else {

               // Reordering successful

               if ( WORK( 2, 1 ).EQ.ZERO ) {

                  // Form C = T22 - lambda*I in WORK(2:N,2:N).

                  DO 20 I = 2, N
                     WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20             CONTINUE
                  N2 = 1
                  NN = N - 1
               } else {

                  // Triangularize the 2 by 2 block by unitary
                  // transformation U = [  cs   i*ss ]
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
               }

               // Estimate norm(inv(C**T))

               EST = ZERO
               KASE = 0
   50          CONTINUE
               CALL SLACN2( NN, WORK( 1, N+2 ), WORK( 1, N+4 ), IWORK, EST, KASE, ISAVE )
               if ( KASE.NE.0 ) {
                  if ( KASE.EQ.1 ) {
                     if ( N2.EQ.1 ) {

                        // Real eigenvalue: solve C**T*x = scale*c.

                        CALL SLAQTR( .TRUE., .TRUE., N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     } else {

                        // Complex eigenvalue: solve
                        // C**T*(p+iq) = scale*(c+id) in real arithmetic.

                        CALL SLAQTR( .TRUE., .FALSE., N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     }
                  } else {
                     if ( N2.EQ.1 ) {

                        // Real eigenvalue: solve C*x = scale*c.

                        CALL SLAQTR( .FALSE., .TRUE., N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )
                     } else {

                        // Complex eigenvalue: solve
                        // C*(p+iq) = scale*(c+id) in real arithmetic.

                        CALL SLAQTR( .FALSE., .FALSE., N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR )

                     }
                  }

                  GO TO 50
               }
            }

            SEP( KS ) = SCALE / MAX( EST, SMLNUM )
            IF( PAIR ) SEP( KS+1 ) = SEP( KS )
         }

         IF( PAIR ) KS = KS + 1

   60 CONTINUE
      RETURN

      // End of STRSNA

      }
