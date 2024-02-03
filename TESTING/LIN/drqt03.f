      SUBROUTINE DRQT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, MINMN, NC;
      double             CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLARNV, DLASET, DORGRQ, DORMRQ
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEED / 1988, 1989, 1990, 1991 /;
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' );
      MINMN = MIN( M, N );

      // Quick return if possible

      if ( MINMN == 0 ) {
         RESULT( 1 ) = ZERO;
         RESULT( 2 ) = ZERO;
         RESULT( 3 ) = ZERO;
         RESULT( 4 ) = ZERO;
         RETURN;
      }

      // Copy the last k rows of the factorization to the array Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if (K > 0 && N > K) CALL DLACPY( 'Full', K, N-K, AF( M-K+1, 1 ), LDA, Q( N-K+1, 1 ), LDA )       IF( K > 1 ) CALL DLACPY( 'Lower', K-1, K-1, AF( M-K+2, N-K+1 ), LDA, Q( N-K+2, N-K+1 ), LDA );

      // Generate the n-by-n matrix Q

      SRNAMT = 'DORGRQ';
      dorgrq(N, N, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK, INFO );

      for (ISIDE = 1; ISIDE <= 2; ISIDE++) { // 30
         if ( ISIDE == 1 ) {
            SIDE = 'L';
            MC = N;
            NC = M;
         } else {
            SIDE = 'R';
            MC = M;
            NC = N;
         }

         // Generate MC by NC matrix C

         for (J = 1; J <= NC; J++) { // 10
            dlarnv(2, ISEED, MC, C( 1, J ) );
         } // 10
         CNORM = DLANGE( '1', MC, NC, C, LDA, RWORK );
         if (CNORM == 0.0) CNORM = ONE;

         for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 20
            if ( ITRANS == 1 ) {
               TRANS = 'N';
            } else {
               TRANS = 'T';
            }

            // Copy C

            dlacpy('Full', MC, NC, C, LDA, CC, LDA );

            // Apply Q or Q' to C

            SRNAMT = 'DORMRQ';
            if (K > 0) CALL DORMRQ( SIDE, TRANS, MC, NC, K, AF( M-K+1, 1 ), LDA, TAU( MINMN-K+1 ), CC, LDA, WORK, LWORK, INFO );

            // Form explicit product and subtract

            if ( LSAME( SIDE, 'L' ) ) {
               dgemm(TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC, LDA );
            } else {
               dgemm('No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC, LDA );
            }

            // Compute error in the difference

            RESID = DLANGE( '1', MC, NC, CC, LDA, RWORK );
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( DBLE( MAX( 1, N ) )*CNORM*EPS );

         } // 20
      } // 30

      RETURN;

      // End of DRQT03

      }
