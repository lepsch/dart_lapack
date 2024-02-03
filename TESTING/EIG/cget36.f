      SUBROUTINE CGET36( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NIN, NINFO;
      REAL               RMAX
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      int                LDT, LWORK;
      PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
      // ..
      // .. Local Scalars ..
      int                I, IFST, ILST, INFO1, INFO2, J, N;
      REAL               EPS, RES
      COMPLEX            CTEMP
      // ..
      // .. Local Arrays ..
      REAL               RESULT( 2 ), RWORK( LDT )
      COMPLEX            DIAG( LDT ), Q( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK )
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHST01, CLACPY, CLASET, CTREXC
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' )
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
      CALL CLACPY( 'F', N, N, TMP, LDT, T1, LDT )
      CALL CLACPY( 'F', N, N, TMP, LDT, T2, LDT )
      RES = ZERO

      // Test without accumulating Q

      CALL CLASET( 'Full', N, N, CZERO, CONE, Q, LDT )
      CALL CTREXC( 'N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1 )
      DO 40 I = 1, N
         DO 30 J = 1, N
            IF( I.EQ.J .AND. Q( I, J ).NE.CONE ) RES = RES + ONE / EPS             IF( I.NE.J .AND. Q( I, J ).NE.CZERO ) RES = RES + ONE / EPS
   30    CONTINUE
   40 CONTINUE

      // Test with accumulating Q

      CALL CLASET( 'Full', N, N, CZERO, CONE, Q, LDT )
      CALL CTREXC( 'V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2 )

      // Compare T1 with T2

      DO 60 I = 1, N
         DO 50 J = 1, N
            IF( T1( I, J ).NE.T2( I, J ) ) RES = RES + ONE / EPS
   50    CONTINUE
   60 CONTINUE
      IF( INFO1.NE.0 .OR. INFO2.NE.0 ) NINFO = NINFO + 1       IF( INFO1.NE.INFO2 ) RES = RES + ONE / EPS

      // Test for successful reordering of T2

      CALL CCOPY( N, TMP, LDT+1, DIAG, 1 )
      IF( IFST.LT.ILST ) THEN
         DO 70 I = IFST + 1, ILST
            CTEMP = DIAG( I )
            DIAG( I ) = DIAG( I-1 )
            DIAG( I-1 ) = CTEMP
   70    CONTINUE
      ELSE IF( IFST.GT.ILST ) THEN
         DO 80 I = IFST - 1, ILST, -1
            CTEMP = DIAG( I+1 )
            DIAG( I+1 ) = DIAG( I )
            DIAG( I ) = CTEMP
   80    CONTINUE
      END IF
      DO 90 I = 1, N
         IF( T2( I, I ).NE.DIAG( I ) ) RES = RES + ONE / EPS
   90 CONTINUE

      // Test for small residual, and orthogonality of Q

      CALL CHST01( N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT )
      RES = RES + RESULT( 1 ) + RESULT( 2 )

      // Test for T2 being in Schur form

      DO 110 J = 1, N - 1
         DO 100 I = J + 1, N
            IF( T2( I, J ).NE.CZERO ) RES = RES + ONE / EPS
  100    CONTINUE
  110 CONTINUE
      IF( RES.GT.RMAX ) THEN
         RMAX = RES
         LMAX = KNT
      END IF
      GO TO 10

      // End of CGET36

      END
