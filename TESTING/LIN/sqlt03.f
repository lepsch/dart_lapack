      SUBROUTINE SQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      REAL               ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, MINMN, NC;
      REAL               CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE;
      // EXTERNAL LSAME, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLARNV, SLASET, SORGQL, SORMQL
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' );
      MINMN = MIN( M, N );

      // Quick return if possible

      if ( MINMN == 0 ) {
         RESULT( 1 ) = ZERO;
         RESULT( 2 ) = ZERO;
         RESULT( 3 ) = ZERO;
         RESULT( 4 ) = ZERO;
         return;
      }

      // Copy the last k columns of the factorization to the array Q

      slaset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      if (K > 0 && M > K) CALL SLACPY( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, M-K+1 ), LDA )       IF( K > 1 ) CALL SLACPY( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, M-K+2 ), LDA );

      // Generate the m-by-m matrix Q

      SRNAMT = 'SORGQL';
      sorgql(M, M, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK, INFO );

      for (ISIDE = 1; ISIDE <= 2; ISIDE++) { // 30
         if ( ISIDE == 1 ) {
            SIDE = 'L';
            MC = M;
            NC = N;
         } else {
            SIDE = 'R';
            MC = N;
            NC = M;
         }

         // Generate MC by NC matrix C

         for (J = 1; J <= NC; J++) { // 10
            slarnv(2, ISEED, MC, C( 1, J ) );
         } // 10
         CNORM = SLANGE( '1', MC, NC, C, LDA, RWORK );
         if (CNORM == 0.0) CNORM = ONE;

         for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 20
            if ( ITRANS == 1 ) {
               TRANS = 'N';
            } else {
               TRANS = 'T';
            }

            // Copy C

            slacpy('Full', MC, NC, C, LDA, CC, LDA );

            // Apply Q or Q' to C

            SRNAMT = 'SORMQL';
            if (K > 0) CALL SORMQL( SIDE, TRANS, MC, NC, K, AF( 1, N-K+1 ), LDA, TAU( MINMN-K+1 ), CC, LDA, WORK, LWORK, INFO );

            // Form explicit product and subtract

            if ( LSAME( SIDE, 'L' ) ) {
               sgemm(TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC, LDA );
            } else {
               sgemm('No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC, LDA );
            }

            // Compute error in the difference

            RESID = SLANGE( '1', MC, NC, CC, LDA, RWORK );
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( REAL( MAX( 1, M ) )*CNORM*EPS );

         } // 20
      } // 30

      return;

      // End of SQLT03

      }
