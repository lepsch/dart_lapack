      void ctrsna(JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      REAL               RWORK( * ), S( * ), SEP( * );
      COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( LDWORK, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0+0 ;
      // ..
      // .. Local Scalars ..
      bool               SOMCON, WANTBH, WANTS, WANTSP;
      String             NORMIN;
      int                I, IERR, IX, J, K, KASE, KS;
      REAL               BIGNUM, EPS, EST, LNRM, RNRM, SCALE, SMLNUM, XNORM;
      COMPLEX            CDUM, PROD;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      COMPLEX            DUMMY( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ICAMAX;
      //- REAL               SCNRM2, SLAMCH;
      //- COMPLEX            CDOTC;
      // EXTERNAL LSAME, ICAMAX, SCNRM2, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLACPY, CLATRS, CSRSCL, CTREXC, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( REAL( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' );
      WANTS = LSAME( JOB, 'E' ) || WANTBH;
      WANTSP = LSAME( JOB, 'V' ) || WANTBH;

      SOMCON = LSAME( HOWMNY, 'S' );

      // Set M to the number of eigenpairs for which condition numbers are
      // to be computed.

      if ( SOMCON ) {
         M = 0;
         for (J = 1; J <= N; J++) { // 10
            if( SELECT( J ) ) M = M + 1;
         } // 10
      } else {
         M = N;
      }

      INFO = 0;
      if ( !WANTS && !WANTSP ) {
         INFO = -1;
      } else if ( !LSAME( HOWMNY, 'A' ) && !SOMCON ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDVL < 1 || ( WANTS && LDVL < N ) ) {
         INFO = -8;
      } else if ( LDVR < 1 || ( WANTS && LDVR < N ) ) {
         INFO = -10;
      } else if ( MM < M ) {
         INFO = -13;
      } else if ( LDWORK < 1 || ( WANTSP && LDWORK < N ) ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('CTRSNA', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( SOMCON ) {
            if( !SELECT( 1 ) ) return;
         }
         if (WANTS) S( 1 ) = ONE;
         IF[WANTSP ) SEP( 1] = ( T( 1, 1 ) ).abs();
         return;
      }

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      KS = 1;
      for (K = 1; K <= N; K++) { // 50

         if ( SOMCON ) {
            if( !SELECT( K ) ) GO TO 50;
         }

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            PROD = CDOTC( N, VR( 1, KS ), 1, VL( 1, KS ), 1 );
            RNRM = SCNRM2( N, VR( 1, KS ), 1 );
            LNRM = SCNRM2( N, VL( 1, KS ), 1 );
            S[KS] = ( PROD ).abs() / ( RNRM*LNRM );

         }

         if ( WANTSP ) {

            // Estimate the reciprocal condition number of the k-th
            // eigenvector.

            // Copy the matrix T to the array WORK and swap the k-th
            // diagonal element to the (1,1) position.

            clacpy('Full', N, N, T, LDT, WORK, LDWORK );
            ctrexc('No Q', N, WORK, LDWORK, DUMMY, 1, K, 1, IERR );

            // Form  C = T22 - lambda*I in WORK(2:N,2:N).

            for (I = 2; I <= N; I++) { // 20
               WORK[I, I] = WORK( I, I ) - WORK( 1, 1 );
            } // 20

            // Estimate a lower bound for the 1-norm of inv(C**H). The 1st
            // and (N+1)th columns of WORK are used to store work vectors.

            SEP[KS] = ZERO;
            EST = ZERO;
            KASE = 0;
            NORMIN = 'N';
            } // 30
            clacn2(N-1, WORK( 1, N+1 ), WORK, EST, KASE, ISAVE );

            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve C**H*x = scale*b

                  clatrs('Upper', 'Conjugate transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR );
               } else {

                  // Solve C*x = scale*b

                  clatrs('Upper', 'No transpose', 'Nonunit', NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, SCALE, RWORK, IERR );
               }
               NORMIN = 'Y';
               if ( SCALE != ONE ) {

                  // Multiply by 1/SCALE if doing so will not cause
                  // overflow.

                  IX = ICAMAX( N-1, WORK, 1 );
                  XNORM = CABS1( WORK( IX, 1 ) );
                  if (SCALE < XNORM*SMLNUM || SCALE == ZERO) GO TO 40;
                  csrscl(N, SCALE, WORK, 1 );
               }
               GO TO 30;
            }

            SEP[KS] = ONE / max( EST, SMLNUM );
         }

         } // 40
         KS = KS + 1;
      } // 50
      return;
      }
