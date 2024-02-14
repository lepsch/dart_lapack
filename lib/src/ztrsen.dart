      void ztrsen(final int JOB, final int COMPQ, final int SELECT, final int N, final Matrix<double> T_, final int LDT, final Matrix<double> Q_, final int LDQ, final int W, final int M, final int S, final int SEP, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final T = T_.dim();
  final Q = Q_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LWORK, M, N;
      double             S, SEP;
      bool               SELECT( * );
      Complex         Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, WANTBH, WANTQ, WANTS, WANTSP;
      int                IERR, K, KASE, KS, LWMIN, N1, N2, NN;
      double             EST, RNORM, SCALE;
      int                ISAVE( 3 );
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             ZLANGE;
      // EXTERNAL lsame, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACN2, ZLACPY, ZTREXC, ZTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Decode and test the input parameters.

      WANTBH = lsame( JOB, 'B' );
      WANTS = lsame( JOB, 'E' ) || WANTBH;
      WANTSP = lsame( JOB, 'V' ) || WANTBH;
      WANTQ = lsame( COMPQ, 'V' );

      // Set M to the number of selected eigenvalues.

      M = 0;
      for (K = 1; K <= N; K++) { // 10
         if( SELECT( K ) ) M = M + 1;
      } // 10

      N1 = M;
      N2 = N - M;
      NN = N1*N2;

      INFO = 0;
      LQUERY = ( LWORK == -1 );

      if ( WANTSP ) {
         LWMIN = max( 1, 2*NN );
      } else if ( lsame( JOB, 'N' ) ) {
         LWMIN = 1;
      } else if ( lsame( JOB, 'E' ) ) {
         LWMIN = max( 1, NN );
      }

      if ( !lsame( JOB, 'N' ) && !WANTS && !WANTSP ) {
         INFO = -1;
      } else if ( !lsame( COMPQ, 'N' ) && !WANTQ ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -8;
      } else if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -14;
      }

      if ( INFO == 0 ) {
         WORK[1] = LWMIN;
      }

      if ( INFO != 0 ) {
         xerbla('ZTRSEN', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == N || M == 0 ) {
         if (WANTS) S = ONE;
         IF( WANTSP ) SEP = ZLANGE( '1', N, N, T, LDT, RWORK );
         GO TO 40;
      }

      // Collect the selected eigenvalues at the top left corner of T.

      KS = 0;
      for (K = 1; K <= N; K++) { // 20
         if ( SELECT( K ) ) {
            KS = KS + 1;

            // Swap the K-th eigenvalue to position KS.

            if (K != KS) ztrexc( COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR );
         }
      } // 20

      if ( WANTS ) {

         // Solve the Sylvester equation for R:

            // T11*R - R*T22 = scale*T12

         zlacpy('F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 );
         ztrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );

         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.

         RNORM = ZLANGE( 'F', N1, N2, WORK, N1, RWORK );
         if ( RNORM == ZERO ) {
            S = ONE;
         } else {
            S = SCALE / ( sqrt( SCALE*SCALE / RNORM+RNORM )* sqrt( RNORM ) );
         }
      }

      if ( WANTSP ) {

         // Estimate sep(T11,T22).

         EST = ZERO;
         KASE = 0;
         } // 30
         zlacn2(NN, WORK( NN+1 ), WORK, EST, KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Solve T11*R - R*T22 = scale*X.

               ztrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            } else {

               // Solve T11**H*R - R*T22**H = scale*X.

               ztrsyl('C', 'C', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            }
            GO TO 30;
         }

         SEP = SCALE / EST;
      }

      } // 40

      // Copy reordered eigenvalues to W.

      for (K = 1; K <= N; K++) { // 50
         W[K] = T( K, K );
      } // 50

      WORK[1] = LWMIN;

      }
