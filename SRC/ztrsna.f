      SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             RWORK( * ), S( * ), SEP( * );
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0+0 ;
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

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) || WANTBH
      WANTSP = LSAME( JOB, 'V' ) || WANTBH

      SOMCON = LSAME( HOWMNY, 'S' )

      // Set M to the number of eigenpairs for which condition numbers are
      // to be computed.

      if ( SOMCON ) {
         M = 0
         for (J = 1; J <= N; J++) { // 10
            IF( SELECT( J ) ) M = M + 1
         } // 10
      } else {
         M = N
      }

      INFO = 0
      if ( !WANTS && !WANTSP ) {
         INFO = -1
      } else if ( !LSAME( HOWMNY, 'A' ) && !SOMCON ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDT < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL < 1 || ( WANTS && LDVL < N ) ) {
         INFO = -8
      } else if ( LDVR < 1 || ( WANTS && LDVR < N ) ) {
         INFO = -10
      } else if ( MM < M ) {
         INFO = -13
      } else if ( LDWORK < 1 || ( WANTSP && LDWORK < N ) ) {
         INFO = -16
      }
      if ( INFO != 0 ) {
         xerbla('ZTRSNA', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( SOMCON ) {
            IF( !SELECT( 1 ) ) RETURN
         }
         if (WANTS) S( 1 ) = ONE          IF( WANTSP ) SEP( 1 ) = ABS( T( 1, 1 ) );
         RETURN
      }

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      KS = 1
      for (K = 1; K <= N; K++) { // 50

         if ( SOMCON ) {
            IF( !SELECT( K ) ) GO TO 50
         }

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            PROD = ZDOTC( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
            RNRM = DZNRM2( N, VR( 1, KS ), 1 )
            LNRM = DZNRM2( N, VL( 1, KS ), 1 )
            S( KS ) = ABS( PROD ) / ( RNRM*LNRM )

         }

         if ( WANTSP ) {

            // Estimate the reciprocal condition number of the k-th
            // eigenvector.

            // Copy the matrix T to the array WORK and swap the k-th
            // diagonal element to the (1,1) position.

            zlacpy('Full', N, N, T, LDT, WORK, LDWORK );
            ztrexc('No Q', N, WORK, LDWORK, DUMMY, 1, K, 1, IERR );

            // Form  C = T22 - lambda*I in WORK(2:N,2:N).

            for (I = 2; I <= N; I++) { // 20
               WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
            } // 20

            // Estimate a lower bound for the 1-norm of inv(C**H). The 1st
            // and (N+1)th columns of WORK are used to store work vectors.

            SEP( KS ) = ZERO
            EST = ZERO
            KASE = 0
            NORMIN = 'N'
            } // 30
            zlacn2(N-1, WORK( 1, N+1 ), WORK, EST, KASE, ISAVE );

            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve C**H*x = scale*b

                  zlatrs('Upper', 'Conjugate transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR );
               } else {

                  // Solve C*x = scale*b

                  zlatrs('Upper', 'No transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR );
               }
               NORMIN = 'Y'
               if ( SCALE != ONE ) {

                  // Multiply by 1/SCALE if doing so will not cause
                  // overflow.

                  IX = IZAMAX( N-1, WORK, 1 )
                  XNORM = CABS1( WORK( IX, 1 ) )
                  if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 40;
                  zdrscl(N, SCALE, WORK, 1 );
               }
               GO TO 30
            }

            SEP( KS ) = ONE / MAX( EST, SMLNUM )
         }

         } // 40
         KS = KS + 1
      } // 50
      RETURN

      // End of ZTRSNA

      }
