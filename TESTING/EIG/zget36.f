      SUBROUTINE ZGET36( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NIN, NINFO;
      double             RMAX;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      int                LDT, LWORK;
      const              LDT = 10, LWORK = 2*LDT*LDT ;
      // ..
      // .. Local Scalars ..
      int                I, IFST, ILST, INFO1, INFO2, J, N;
      double             EPS, RES;
      COMPLEX*16         CTEMP
      // ..
      // .. Local Arrays ..
      double             RESULT( 2 ), RWORK( LDT );
      COMPLEX*16         DIAG( LDT ), Q( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK )
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZHST01, ZLACPY, ZLASET, ZTREXC
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'P' )
      RMAX = ZERO
      LMAX = 0
      KNT = 0
      NINFO = 0

      // Read input data until N=0

   10 CONTINUE
      READ( NIN, FMT = * )N, IFST, ILST
      IF( N.EQ.0 ) RETURN
      KNT = KNT + 1
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
      zlacpy('F', N, N, TMP, LDT, T1, LDT );
      zlacpy('F', N, N, TMP, LDT, T2, LDT );
      RES = ZERO

      // Test without accumulating Q

      zlaset('Full', N, N, CZERO, CONE, Q, LDT );
      ztrexc('N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1 );
      DO 40 I = 1, N
         DO 30 J = 1, N
            IF( I.EQ.J .AND. Q( I, J ).NE.CONE ) RES = RES + ONE / EPS             IF( I.NE.J .AND. Q( I, J ).NE.CZERO ) RES = RES + ONE / EPS
   30    CONTINUE
   40 CONTINUE

      // Test with accumulating Q

      zlaset('Full', N, N, CZERO, CONE, Q, LDT );
      ztrexc('V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2 );

      // Compare T1 with T2

      DO 60 I = 1, N
         DO 50 J = 1, N
            IF( T1( I, J ).NE.T2( I, J ) ) RES = RES + ONE / EPS
   50    CONTINUE
   60 CONTINUE
      IF( INFO1.NE.0 .OR. INFO2.NE.0 ) NINFO = NINFO + 1       IF( INFO1.NE.INFO2 ) RES = RES + ONE / EPS

      // Test for successful reordering of T2

      zcopy(N, TMP, LDT+1, DIAG, 1 );
      if ( IFST.LT.ILST ) {
         DO 70 I = IFST + 1, ILST
            CTEMP = DIAG( I )
            DIAG( I ) = DIAG( I-1 )
            DIAG( I-1 ) = CTEMP
   70    CONTINUE
      } else if ( IFST.GT.ILST ) {
         DO 80 I = IFST - 1, ILST, -1
            CTEMP = DIAG( I+1 )
            DIAG( I+1 ) = DIAG( I )
            DIAG( I ) = CTEMP
   80    CONTINUE
      }
      DO 90 I = 1, N
         IF( T2( I, I ).NE.DIAG( I ) ) RES = RES + ONE / EPS
   90 CONTINUE

      // Test for small residual, and orthogonality of Q

      zhst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT );
      RES = RES + RESULT( 1 ) + RESULT( 2 )

      // Test for T2 being in Schur form

      DO 110 J = 1, N - 1
         DO 100 I = J + 1, N
            IF( T2( I, J ).NE.CZERO ) RES = RES + ONE / EPS
  100    CONTINUE
  110 CONTINUE
      if ( RES.GT.RMAX ) {
         RMAX = RES
         LMAX = KNT
      }
      GO TO 10

      // End of ZGET36

      }
