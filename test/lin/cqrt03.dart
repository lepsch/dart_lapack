      void cqrt03(M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * );
      Complex            AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      REAL               CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLARNV, CLASET, CUNGQR, CUNMQR
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
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

      // Copy the first k columns of the factorization to the array Q

      claset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      clacpy('Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the m-by-m matrix Q

      SRNAMT = 'CUNGQR';
      cungqr(M, M, K, Q, LDA, TAU, WORK, LWORK, INFO );

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

            SRNAMT = 'CUNMQR';
            cunmqr(SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO );

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

      return;
      }
