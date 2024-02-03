      SUBROUTINE CTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S, SEP, WORK, LWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LWORK, M, N;
      REAL               S, SEP;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTBH, WANTQ, WANTS, WANTSP;
      int                IERR, K, KASE, KS, LWMIN, N1, N2, NN;
      REAL               EST, RNORM, SCALE;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      REAL               RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SROUNDUP_LWORK;
      // EXTERNAL LSAME, CLANGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLACPY, CTREXC, CTRSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      WANTBH = LSAME( JOB, 'B' );
      WANTS = LSAME( JOB, 'E' ) || WANTBH;
      WANTSP = LSAME( JOB, 'V' ) || WANTBH;
      WANTQ = LSAME( COMPQ, 'V' );

      // Set M to the number of selected eigenvalues.

      M = 0;
      for (K = 1; K <= N; K++) { // 10
         IF( SELECT( K ) ) M = M + 1;
      } // 10

      N1 = M;
      N2 = N - M;
      NN = N1*N2;

      INFO = 0;
      LQUERY = ( LWORK == -1 );

      if ( WANTSP ) {
         LWMIN = MAX( 1, 2*NN );
      } else if ( LSAME( JOB, 'N' ) ) {
         LWMIN = 1;
      } else if ( LSAME( JOB, 'E' ) ) {
         LWMIN = MAX( 1, NN );
      }

      if ( !LSAME( JOB, 'N' ) && !WANTS && !WANTSP ) {
         INFO = -1;
      } else if ( !LSAME( COMPQ, 'N' ) && !WANTQ ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < MAX( 1, N ) ) {
         INFO = -6;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -8;
      } else if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -14;
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN);
      }

      if ( INFO != 0 ) {
         xerbla('CTRSEN', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == N || M == 0 ) {
         if (WANTS) S = ONE          IF( WANTSP ) SEP = CLANGE( '1', N, N, T, LDT, RWORK );
         GO TO 40;
      }

      // Collect the selected eigenvalues at the top left corner of T.

      KS = 0;
      for (K = 1; K <= N; K++) { // 20
         if ( SELECT( K ) ) {
            KS = KS + 1;

            // Swap the K-th eigenvalue to position KS.

            if (K != KS) CALL CTREXC( COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR );
         }
      } // 20

      if ( WANTS ) {

         // Solve the Sylvester equation for R:

            // T11*R - R*T22 = scale*T12

         clacpy('F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 );
         ctrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );

         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.

         RNORM = CLANGE( 'F', N1, N2, WORK, N1, RWORK );
         if ( RNORM == ZERO ) {
            S = ONE;
         } else {
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* SQRT( RNORM ) );
         }
      }

      if ( WANTSP ) {

         // Estimate sep(T11,T22).

         EST = ZERO;
         KASE = 0;
         } // 30
         clacn2(NN, WORK( NN+1 ), WORK, EST, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Solve T11*R - R*T22 = scale*X.

               ctrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            } else {

               // Solve T11**H*R - R*T22**H = scale*X.

               ctrsyl('C', 'C', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            }
            GO TO 30;
         }

         SEP = SCALE / EST;
      }

      } // 40

      // Copy reordered eigenvalues to W.

      for (K = 1; K <= N; K++) { // 50
         W( K ) = T( K, K );
      } // 50

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN);

      return;

      // End of CTRSEN

      }
