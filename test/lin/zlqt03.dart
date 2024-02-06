      void zlqt03(M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double             RESULT( * ), RWORK( * );
      Complex         AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      double             CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZLARNV, ZLASET, ZUNGLQ, ZUNMLQ
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
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

      zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      zlacpy('Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

     srnamc.SRNAMT = 'ZUNGLQ';
      zunglq(N, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

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
            zlarnv(2, ISEED, MC, C( 1, J ) );
         } // 10
         CNORM = ZLANGE( '1', MC, NC, C, LDA, RWORK );
         if (CNORM == ZERO) CNORM = ONE;

         for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 20
            if ( ITRANS == 1 ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            // Copy C

            zlacpy('Full', MC, NC, C, LDA, CC, LDA );

            // Apply Q or Q' to C

           srnamc.SRNAMT = 'ZUNMLQ';
            zunmlq(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO );

            // Form explicit product and subtract

            if ( lsame( SIDE, 'L' ) ) {
               zgemm(TRANS, 'No transpose', MC, NC, MC, DCMPLX( -ONE ), Q, LDA, C, LDA, DCMPLX( ONE ), CC, LDA );
            } else {
               zgemm('No transpose', TRANS, MC, NC, NC, DCMPLX( -ONE ), C, LDA, Q, LDA, DCMPLX( ONE ), CC, LDA );
            }

            // Compute error in the difference

            RESID = ZLANGE( '1', MC, NC, CC, LDA, RWORK );
            RESULT[( ISIDE-1 )*2+ITRANS] = RESID / ( (max( 1, N )).toDouble()*CNORM*EPS );

         } // 20
      } // 30

      return;
      }
