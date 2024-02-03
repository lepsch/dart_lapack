      SUBROUTINE DGET40( RMAX, LMAX, NINFO, KNT, NIN )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                KNT, LMAX, NIN;
      double             RMAX;
*     ..
*     .. Array Arguments ..
      int                NINFO( 2 );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      int                LDT, LWORK;
      PARAMETER          ( LDT = 10, LWORK = 100 + 4*LDT + 16 )
*     ..
*     .. Local Scalars ..
      int                I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, ILST2, ILSTSV, J, LOC, N;
      double             EPS, RES;
*     ..
*     .. Local Arrays ..
      double             Q( LDT, LDT ), Z( LDT, LDT ), RESULT( 4 ), T( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), S( LDT, LDT ), S1( LDT, LDT ), S2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK );
*     ..
*     .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
*     ..
*     .. External Subroutines ..
      // EXTERNAL DHST01, DLACPY, DLASET, DTGEXC
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'P' )
      RMAX = ZERO
      LMAX = 0
      KNT = 0
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
*
*     Read input data until N=0
*
   10 CONTINUE
      READ( NIN, FMT = * )N, IFST, ILST
      IF( N.EQ.0 ) RETURN
      KNT = KNT + 1
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
      CALL DLACPY( 'F', N, N, TMP, LDT, T, LDT )
      CALL DLACPY( 'F', N, N, TMP, LDT, T1, LDT )
      CALL DLACPY( 'F', N, N, TMP, LDT, T2, LDT )
      DO 25 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   25 CONTINUE
      CALL DLACPY( 'F', N, N, TMP, LDT, S, LDT )
      CALL DLACPY( 'F', N, N, TMP, LDT, S1, LDT )
      CALL DLACPY( 'F', N, N, TMP, LDT, S2, LDT )
      IFSTSV = IFST
      ILSTSV = ILST
      IFST1 = IFST
      ILST1 = ILST
      IFST2 = IFST
      ILST2 = ILST
      RES = ZERO
*
*     Test without accumulating Q and Z
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
      CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDT )
      CALL DTGEXC( .FALSE., .FALSE., N, T1, LDT, S1, LDT, Q, LDT, Z, LDT, IFST1, ILST1, WORK, LWORK, NINFO ( 1 ) )
      DO 40 I = 1, N
         DO 30 J = 1, N
            IF( I.EQ.J .AND. Q( I, J ).NE.ONE ) RES = RES + ONE / EPS             IF( I.NE.J .AND. Q( I, J ).NE.ZERO ) RES = RES + ONE / EPS             IF( I.EQ.J .AND. Z( I, J ).NE.ONE ) RES = RES + ONE / EPS             IF( I.NE.J .AND. Z( I, J ).NE.ZERO ) RES = RES + ONE / EPS
   30    CONTINUE
   40 CONTINUE
*
*     Test with accumulating Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDT )
      CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDT )
      CALL DTGEXC( .TRUE., .TRUE., N, T2, LDT, S2, LDT, Q, LDT, Z, LDT, IFST2, ILST2, WORK, LWORK, NINFO ( 2 ) )
*
*     Compare T1 with T2 and S1 with S2
*
      DO 60 I = 1, N
         DO 50 J = 1, N
            IF( T1( I, J ).NE.T2( I, J ) ) RES = RES + ONE / EPS             IF( S1( I, J ).NE.S2( I, J ) ) RES = RES + ONE / EPS
   50    CONTINUE
   60 CONTINUE
      IF( IFST1.NE.IFST2 ) RES = RES + ONE / EPS       IF( ILST1.NE.ILST2 ) RES = RES + ONE / EPS       IF( NINFO( 1 ).NE.NINFO( 2 ) ) RES = RES + ONE / EPS
*
*     Test orthogonality of Q and Z and backward error on T2 and S2
*
      CALL DGET51( 1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 1 ) )       CALL DGET51( 1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 2 ) )       CALL DGET51( 3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK, RESULT( 3 ) )       CALL DGET51( 3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK, RESULT( 4 ) )
      RES = RES + RESULT( 1 ) + RESULT( 2 ) + RESULT( 3 ) + RESULT( 4 )
*
*     Read next matrix pair
*
      GO TO 10
*
*     End of DGET40
*
      END
