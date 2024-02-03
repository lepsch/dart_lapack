      SUBROUTINE DORT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             U( LDU, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN;
      double             EPS, TMP;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT, DLAMCH, DLANSY;
      // EXTERNAL LSAME, DDOT, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESID = ZERO

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      EPS = DLAMCH( 'Precision' )
      if ( M.LT.N .OR. ( M.EQ.N .AND. LSAME( ROWCOL, 'R' ) ) ) {
         TRANSU = 'N'
         K = N
      } else {
         TRANSU = 'T'
         K = M
      }
      MNMIN = MIN( M, N )

      if ( ( MNMIN+1 )*MNMIN.LE.LWORK ) {
         LDWORK = MNMIN
      } else {
         LDWORK = 0
      }
      if ( LDWORK.GT.0 ) {

         // Compute I - U*U' or I - U'*U.

         dlaset('Upper', MNMIN, MNMIN, ZERO, ONE, WORK, LDWORK );
         dsyrk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK );

         // Compute norm( I - U*U' ) / ( K * EPS ) .

         RESID = DLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, WORK( LDWORK*MNMIN+1 ) )
         RESID = ( RESID / DBLE( K ) ) / EPS
      } else if ( TRANSU.EQ.'T' ) {

         // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

         DO 20 J = 1, N
            DO 10 I = 1, J
               if ( I.NE.J ) {
                  TMP = ZERO
               } else {
                  TMP = ONE
               }
               TMP = TMP - DDOT( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, ABS( TMP ) )
   10       CONTINUE
   20    CONTINUE
         RESID = ( RESID / DBLE( M ) ) / EPS
      } else {

         // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

         DO 40 J = 1, M
            DO 30 I = 1, J
               if ( I.NE.J ) {
                  TMP = ZERO
               } else {
                  TMP = ONE
               }
               TMP = TMP - DDOT( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, ABS( TMP ) )
   30       CONTINUE
   40    CONTINUE
         RESID = ( RESID / DBLE( N ) ) / EPS
      }
      RETURN

      // End of DORT01

      }
