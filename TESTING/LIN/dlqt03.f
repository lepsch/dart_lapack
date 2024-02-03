      SUBROUTINE DLQT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT )

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
      double             ONE;
      const              ONE = 1.0D0 ;
      double             ROGUE;
      const              ROGUE = -1.0D+10 ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      double             CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLARNV, DLASET, DORGLQ, DORMLQ
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
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

      EPS = DLAMCH( 'Epsilon' )

      // Copy the first k rows of the factorization to the array Q

      CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      CALL DLACPY( 'Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )

      // Generate the n-by-n matrix Q

      SRNAMT = 'DORGLQ'
      CALL DORGLQ( N, N, K, Q, LDA, TAU, WORK, LWORK, INFO )

      DO 30 ISIDE = 1, 2
         if ( ISIDE.EQ.1 ) {
            SIDE = 'L'
            MC = N
            NC = M
         } else {
            SIDE = 'R'
            MC = M
            NC = N
         }

         // Generate MC by NC matrix C

         DO 10 J = 1, NC
            CALL DLARNV( 2, ISEED, MC, C( 1, J ) )
   10    CONTINUE
         CNORM = DLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.0.0D0 ) CNORM = ONE

         DO 20 ITRANS = 1, 2
            if ( ITRANS.EQ.1 ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            // Copy C

            CALL DLACPY( 'Full', MC, NC, C, LDA, CC, LDA )

            // Apply Q or Q' to C

            SRNAMT = 'DORMLQ'
            CALL DORMLQ( SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO )

            // Form explicit product and subtract

            if ( LSAME( SIDE, 'L' ) ) {
               CALL DGEMM( TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC, LDA )
            } else {
               CALL DGEMM( 'No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC, LDA )
            }

            // Compute error in the difference

            RESID = DLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( DBLE( MAX( 1, N ) )*CNORM*EPS )

   20    CONTINUE
   30 CONTINUE

      RETURN

      // End of DLQT03

      }
