      SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S, SEP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, JOB;
      int                INFO, LDQ, LDT, LWORK, M, N;
      double             S, SEP;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      COMPLEX*16         Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
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

      // Decode and test the input parameters.

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )

      // Set M to the number of selected eigenvalues.

      M = 0
      for (K = 1; K <= N; K++) { // 10
         IF( SELECT( K ) ) M = M + 1
      } // 10

      N1 = M
      N2 = N - M
      NN = N1*N2

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )

      if ( WANTSP ) {
         LWMIN = MAX( 1, 2*NN )
      } else if ( LSAME( JOB, 'N' ) ) {
         LWMIN = 1
      } else if ( LSAME( JOB, 'E' ) ) {
         LWMIN = MAX( 1, NN )
      }

      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP ) {
         INFO = -1
      } else if ( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -8
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -14
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 ) = LWMIN
      }

      if ( INFO.NE.0 ) {
         xerbla('ZTRSEN', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.N .OR. M.EQ.0 ) {
         if (WANTS) S = ONE          IF( WANTSP ) SEP = ZLANGE( '1', N, N, T, LDT, RWORK );
         GO TO 40
      }

      // Collect the selected eigenvalues at the top left corner of T.

      KS = 0
      for (K = 1; K <= N; K++) { // 20
         if ( SELECT( K ) ) {
            KS = KS + 1

            // Swap the K-th eigenvalue to position KS.

            if (K.NE.KS) CALL ZTREXC( COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR );
         }
      } // 20

      if ( WANTS ) {

         // Solve the Sylvester equation for R:

            // T11*R - R*T22 = scale*T12

         zlacpy('F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 );
         ztrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );

         // Estimate the reciprocal of the condition number of the cluster
         // of eigenvalues.

         RNORM = ZLANGE( 'F', N1, N2, WORK, N1, RWORK )
         if ( RNORM.EQ.ZERO ) {
            S = ONE
         } else {
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* SQRT( RNORM ) )
         }
      }

      if ( WANTSP ) {

         // Estimate sep(T11,T22).

         EST = ZERO
         KASE = 0
         } // 30
         zlacn2(NN, WORK( NN+1 ), WORK, EST, KASE, ISAVE );
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.1 ) {

               // Solve T11*R - R*T22 = scale*X.

               ztrsyl('N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            } else {

               // Solve T11**H*R - R*T22**H = scale*X.

               ztrsyl('C', 'C', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, IERR );
            }
            GO TO 30
         }

         SEP = SCALE / EST
      }

      } // 40

      // Copy reordered eigenvalues to W.

      for (K = 1; K <= N; K++) { // 50
         W( K ) = T( K, K )
      } // 50

      WORK( 1 ) = LWMIN

      RETURN

      // End of ZTRSEN

      }
