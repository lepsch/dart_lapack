      SUBROUTINE ZQRT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      COMPLEX*16         AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         ROGUE
      const              ROGUE = ( -1.0D+10, -1.0D+10 ) ;
      // ..
      // .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      double             CNORM, EPS, RESID;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZLARNV, ZLASET, ZUNGQR, ZUNMQR
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
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

      // Copy the first k columns of the factorization to the array Q

      CALL ZLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      CALL ZLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )

      // Generate the m-by-m matrix Q

      SRNAMT = 'ZUNGQR'
      CALL ZUNGQR( M, M, K, Q, LDA, TAU, WORK, LWORK, INFO )

      DO 30 ISIDE = 1, 2
         IF( ISIDE.EQ.1 ) THEN
            SIDE = 'L'
            MC = M
            NC = N
         } else {
            SIDE = 'R'
            MC = N
            NC = M
         END IF

         // Generate MC by NC matrix C

         DO 10 J = 1, NC
            CALL ZLARNV( 2, ISEED, MC, C( 1, J ) )
   10    CONTINUE
         CNORM = ZLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.ZERO ) CNORM = ONE

         DO 20 ITRANS = 1, 2
            IF( ITRANS.EQ.1 ) THEN
               TRANS = 'N'
            } else {
               TRANS = 'C'
            END IF

            // Copy C

            CALL ZLACPY( 'Full', MC, NC, C, LDA, CC, LDA )

            // Apply Q or Q' to C

            SRNAMT = 'ZUNMQR'
            CALL ZUNMQR( SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO )

            // Form explicit product and subtract

            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL ZGEMM( TRANS, 'No transpose', MC, NC, MC, DCMPLX( -ONE ), Q, LDA, C, LDA, DCMPLX( ONE ), CC, LDA )
            } else {
               CALL ZGEMM( 'No transpose', TRANS, MC, NC, NC, DCMPLX( -ONE ), C, LDA, Q, LDA, DCMPLX( ONE ), CC, LDA )
            END IF

            // Compute error in the difference

            RESID = ZLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( DBLE( MAX( 1, M ) )*CNORM*EPS )

   20    CONTINUE
   30 CONTINUE

      RETURN

      // End of ZQRT03

      }
