      SUBROUTINE SGET36( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NIN;
      REAL               RMAX
      // ..
      // .. Array Arguments ..
      int                NINFO( 3 );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      int                LDT, LWORK;
      PARAMETER          ( LDT = 10, LWORK = 2*LDT*LDT )
      // ..
      // .. Local Scalars ..
      int                I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, ILST2, ILSTSV, INFO1, INFO2, J, LOC, N;
      REAL               EPS, RES
      // ..
      // .. Local Arrays ..
      REAL               Q( LDT, LDT ), RESULT( 2 ), T1( LDT, LDT ), T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK )
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SHST01, SLACPY, SLASET, STREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' )
      RMAX = ZERO
      LMAX = 0
      KNT = 0
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
      NINFO( 3 ) = 0

      // Read input data until N=0

   10 CONTINUE
      READ( NIN, FMT = * )N, IFST, ILST
      IF( N.EQ.0 ) RETURN
      KNT = KNT + 1
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
      CALL SLACPY( 'F', N, N, TMP, LDT, T1, LDT )
      CALL SLACPY( 'F', N, N, TMP, LDT, T2, LDT )
      IFSTSV = IFST
      ILSTSV = ILST
      IFST1 = IFST
      ILST1 = ILST
      IFST2 = IFST
      ILST2 = ILST
      RES = ZERO

      // Test without accumulating Q

      CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
      CALL STREXC( 'N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1 )
      DO 40 I = 1, N
         DO 30 J = 1, N
            IF( I.EQ.J .AND. Q( I, J ).NE.ONE ) RES = RES + ONE / EPS             IF( I.NE.J .AND. Q( I, J ).NE.ZERO ) RES = RES + ONE / EPS
   30    CONTINUE
   40 CONTINUE

      // Test with accumulating Q

      CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
      CALL STREXC( 'V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2 )

      // Compare T1 with T2

      DO 60 I = 1, N
         DO 50 J = 1, N
            IF( T1( I, J ).NE.T2( I, J ) ) RES = RES + ONE / EPS
   50    CONTINUE
   60 CONTINUE
      IF( IFST1.NE.IFST2 ) RES = RES + ONE / EPS       IF( ILST1.NE.ILST2 ) RES = RES + ONE / EPS       IF( INFO1.NE.INFO2 ) RES = RES + ONE / EPS

      // Test for successful reordering of T2

      IF( INFO2.NE.0 ) THEN
         NINFO( INFO2 ) = NINFO( INFO2 ) + 1
      ELSE
         IF( ABS( IFST2-IFSTSV ).GT.1 ) RES = RES + ONE / EPS          IF( ABS( ILST2-ILSTSV ).GT.1 ) RES = RES + ONE / EPS
      END IF

      // Test for small residual, and orthogonality of Q

      CALL SHST01( N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RESULT )
      RES = RES + RESULT( 1 ) + RESULT( 2 )

      // Test for T2 being in Schur form

      LOC = 1
   70 CONTINUE
      IF( T2( LOC+1, LOC ).NE.ZERO ) THEN

         // 2 by 2 block

         IF( T2( LOC, LOC+1 ).EQ.ZERO .OR. T2( LOC, LOC ).NE. T2( LOC+1, LOC+1 ) .OR. SIGN( ONE, T2( LOC, LOC+1 ) ).EQ. SIGN( ONE, T2( LOC+1, LOC ) ) )RES = RES + ONE / EPS
         DO 80 I = LOC + 2, N
            IF( T2( I, LOC ).NE.ZERO ) RES = RES + ONE / RES             IF( T2( I, LOC+1 ).NE.ZERO ) RES = RES + ONE / RES
   80    CONTINUE
         LOC = LOC + 2
      ELSE

         // 1 by 1 block

         DO 90 I = LOC + 1, N
            IF( T2( I, LOC ).NE.ZERO ) RES = RES + ONE / RES
   90    CONTINUE
         LOC = LOC + 1
      END IF
      IF( LOC.LT.N ) GO TO 70
      IF( RES.GT.RMAX ) THEN
         RMAX = RES
         LMAX = KNT
      END IF
      GO TO 10

      // End of SGET36

      END
