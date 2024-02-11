      void dlqt03(final int M, final int N, final int K, final int AF, final int C, final int CC, final int Q, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double             AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      double             CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLARNV, DLASET, DORGLQ, DORMLQ
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = dlamch( 'Epsilon' );

      // Copy the first k rows of the factorization to the array Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      dlacpy('Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

     srnamc.SRNAMT = 'DORGLQ';
      dorglq(N, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

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
         CNORM = dlange( '1', MC, NC, C, LDA, RWORK );
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

           srnamc.SRNAMT = 'DORMLQ';
            dormlq(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO );

            // Form explicit product and subtract

            if ( lsame( SIDE, 'L' ) ) {
               dgemm(TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC, LDA );
            } else {
               dgemm('No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC, LDA );
            }

            // Compute error in the difference

            RESID = dlange( '1', MC, NC, CC, LDA, RWORK );
            RESULT[( ISIDE-1 )*2+ITRANS] = RESID / ( (max( 1, N )).toDouble()*CNORM*EPS );

         } // 20
      } // 30

      }
