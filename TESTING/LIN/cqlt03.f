      SUBROUTINE CQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            ROGUE
      const              ROGUE = ( -1.0E+10, -1.0E+10 ) ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, MINMN, NC;
      REAL               CNORM, EPS, RESID
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLARNV, CLASET, CUNGQL, CUNMQL
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEED / 1988, 1989, 1990, 1991 /
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )
      MINMN = MIN( M, N )

      // Quick return if possible

      if ( MINMN.EQ.0 ) {
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RESULT( 3 ) = ZERO
         RESULT( 4 ) = ZERO
         RETURN
      ENDIF

      // Copy the last k columns of the factorization to the array Q

      claset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      IF( K.GT.0 .AND. M.GT.K ) CALL CLACPY( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, M-K+1 ), LDA )       IF( K.GT.1 ) CALL CLACPY( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, M-K+2 ), LDA )

      // Generate the m-by-m matrix Q

      SRNAMT = 'CUNGQL'
      cungql(M, M, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK, INFO );

      for (ISIDE = 1; ISIDE <= 2; ISIDE++) { // 30
         if ( ISIDE.EQ.1 ) {
            SIDE = 'L'
            MC = M
            NC = N
         } else {
            SIDE = 'R'
            MC = N
            NC = M
         }

         // Generate MC by NC matrix C

         for (J = 1; J <= NC; J++) { // 10
            clarnv(2, ISEED, MC, C( 1, J ) );
         } // 10
         CNORM = CLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.ZERO ) CNORM = ONE

         for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 20
            if ( ITRANS.EQ.1 ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            // Copy C

            clacpy('Full', MC, NC, C, LDA, CC, LDA );

            // Apply Q or Q' to C

            SRNAMT = 'CUNMQL'
            IF( K.GT.0 ) CALL CUNMQL( SIDE, TRANS, MC, NC, K, AF( 1, N-K+1 ), LDA, TAU( MINMN-K+1 ), CC, LDA, WORK, LWORK, INFO )

            // Form explicit product and subtract

            if ( LSAME( SIDE, 'L' ) ) {
               cgemm(TRANS, 'No transpose', MC, NC, MC, CMPLX( -ONE ), Q, LDA, C, LDA, CMPLX( ONE ), CC, LDA );
            } else {
               cgemm('No transpose', TRANS, MC, NC, NC, CMPLX( -ONE ), C, LDA, Q, LDA, CMPLX( ONE ), CC, LDA );
            }

            // Compute error in the difference

            RESID = CLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( REAL( MAX( 1, M ) )*CNORM*EPS )

         } // 20
      } // 30

      RETURN

      // End of CQLT03

      }
