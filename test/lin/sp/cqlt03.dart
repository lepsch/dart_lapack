      void cqlt03(M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               RESULT( * ), RWORK( * );
      Complex            AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, MINMN, NC;
      double               CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLARNV, CLASET, CUNGQL, CUNMQL
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = SLAMCH( 'Epsilon' );
      MINMN = min( M, N );

      // Quick return if possible

      if ( MINMN == 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         RESULT[3] = ZERO;
         RESULT[4] = ZERO;
         return;
      }

      // Copy the last k columns of the factorization to the array Q

      claset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      if (K > 0 && M > K) clacpy( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, M-K+1 ), LDA );
      IF( K > 1 ) clacpy( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, M-K+2 ), LDA );

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'CUNGQL';
      cungql(M, M, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK, INFO );

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
            clarnv(2, ISEED, MC, C( 1, J ) );
         } // 10
         CNORM = CLANGE( '1', MC, NC, C, LDA, RWORK );
         if (CNORM == ZERO) CNORM = ONE;

         for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 20
            if ( ITRANS == 1 ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            // Copy C

            clacpy('Full', MC, NC, C, LDA, CC, LDA );

            // Apply Q or Q' to C

           srnamc.SRNAMT = 'CUNMQL';
            if (K > 0) cunmql( SIDE, TRANS, MC, NC, K, AF( 1, N-K+1 ), LDA, TAU( MINMN-K+1 ), CC, LDA, WORK, LWORK, INFO );

            // Form explicit product and subtract

            if ( lsame( SIDE, 'L' ) ) {
               cgemm(TRANS, 'No transpose', MC, NC, MC, CMPLX( -ONE ), Q, LDA, C, LDA, CMPLX( ONE ), CC, LDA );
            } else {
               cgemm('No transpose', TRANS, MC, NC, NC, CMPLX( -ONE ), C, LDA, Q, LDA, CMPLX( ONE ), CC, LDA );
            }

            // Compute error in the difference

            RESID = CLANGE( '1', MC, NC, CC, LDA, RWORK );
            RESULT[( ISIDE-1 )*2+ITRANS] = RESID / ( REAL( max( 1, M ) )*CNORM*EPS );

         } // 20
      } // 30

      }
