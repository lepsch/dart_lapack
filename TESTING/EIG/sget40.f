      SUBROUTINE SGET40( RMAX, LMAX, NINFO, KNT, NIN );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NIN;
      REAL               RMAX;
      // ..
      // .. Array Arguments ..
      int                NINFO( 2 );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                LDT, LWORK;
      const              LDT = 10, LWORK = 100 + 4*LDT + 16 ;
      // ..
      // .. Local Scalars ..
      int                I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, ILST2, ILSTSV, J, LOC, N;
      REAL               EPS, RES;
      // ..
      // .. Local Arrays ..
      REAL               Q( LDT, LDT ), Z( LDT, LDT ), RESULT( 4 ), T( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), S( LDT, LDT ), S1( LDT, LDT ), S2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK );
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGET51, SLACPY, SLASET, STGEXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' );
      RMAX = ZERO;
      LMAX = 0;
      KNT = 0;
      NINFO( 1 ) = 0;
      NINFO( 2 ) = 0;

      // Read input data until N=0

      } // 10
      READ( NIN, FMT = * )N, IFST, ILST;
      if (N == 0) RETURN;
      KNT = KNT + 1;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 20
      slacpy('F', N, N, TMP, LDT, T, LDT );
      slacpy('F', N, N, TMP, LDT, T1, LDT );
      slacpy('F', N, N, TMP, LDT, T2, LDT );
      for (I = 1; I <= N; I++) { // 25
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 25
      slacpy('F', N, N, TMP, LDT, S, LDT );
      slacpy('F', N, N, TMP, LDT, S1, LDT );
      slacpy('F', N, N, TMP, LDT, S2, LDT );
      IFSTSV = IFST;
      ILSTSV = ILST;
      IFST1 = IFST;
      ILST1 = ILST;
      IFST2 = IFST;
      ILST2 = ILST;
      RES = ZERO;

      // Test without accumulating Q and Z

      slaset('Full', N, N, ZERO, ONE, Q, LDT );
      slaset('Full', N, N, ZERO, ONE, Z, LDT );
      stgexc( false , false , N, T1, LDT, S1, LDT, Q, LDT, Z, LDT, IFST1, ILST1, WORK, LWORK, NINFO( 1 ) );
      for (I = 1; I <= N; I++) { // 40
         for (J = 1; J <= N; J++) { // 30
            IF( I == J && Q( I, J ) != ONE ) RES = RES + ONE / EPS             IF( I != J && Q( I, J ) != ZERO ) RES = RES + ONE / EPS             IF( I == J && Z( I, J ) != ONE ) RES = RES + ONE / EPS             IF( I != J && Z( I, J ) != ZERO ) RES = RES + ONE / EPS;
         } // 30
      } // 40

      // Test with accumulating Q

      slaset('Full', N, N, ZERO, ONE, Q, LDT );
      slaset('Full', N, N, ZERO, ONE, Z, LDT );
      stgexc( true , true , N, T2, LDT, S2, LDT, Q, LDT, Z, LDT, IFST2, ILST2, WORK, LWORK, NINFO( 2 ) );

      // Compare T1 with T2 and S1 with S2

      for (I = 1; I <= N; I++) { // 60
         for (J = 1; J <= N; J++) { // 50
            IF( T1( I, J ) != T2( I, J ) ) RES = RES + ONE / EPS             IF( S1( I, J ) != S2( I, J ) ) RES = RES + ONE / EPS;
         } // 50
      } // 60
      if (IFST1 != IFST2) RES = RES + ONE / EPS       IF( ILST1 != ILST2 ) RES = RES + ONE / EPS       IF( NINFO( 1 ) != NINFO( 2 ) ) RES = RES + ONE / EPS;

      // Test orthogonality of Q and Z and backward error on T2 and S2

      sget51(1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 1 ) )       CALL SGET51( 1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK, RESULT( 2 ) )       CALL SGET51( 3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK, RESULT( 3 ) )       CALL SGET51( 3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK, RESULT( 4 ) );
      RES = RES + RESULT( 1 ) + RESULT( 2 ) + RESULT( 3 ) + RESULT( 4 );

      // Read next matrix pair

      GO TO 10;

      // End of SGET40

      }
