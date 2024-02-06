      void sget36(RMAX, LMAX, NINFO, KNT, NIN ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KNT, LMAX, NIN;
      double               RMAX;
      int                NINFO( 3 );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                LDT, LWORK;
      const              LDT = 10, LWORK = 2*LDT*LDT ;
      int                I, IFST, IFST1, IFST2, IFSTSV, ILST, ILST1, ILST2, ILSTSV, INFO1, INFO2, J, LOC, N;
      double               EPS, RES;
      double               Q( LDT, LDT ), RESULT( 2 ), T1( LDT, LDT ), T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SHST01, SLACPY, SLASET, STREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN

      EPS = SLAMCH( 'P' );
      RMAX = ZERO;
      LMAX = 0;
      KNT = 0;
      NINFO[1] = 0;
      NINFO[2] = 0;
      NINFO[3] = 0;

      // Read input data until N=0

      } // 10
      READ( NIN, FMT = * )N, IFST, ILST;
      if (N == 0) return;
      KNT = KNT + 1;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 20
      slacpy('F', N, N, TMP, LDT, T1, LDT );
      slacpy('F', N, N, TMP, LDT, T2, LDT );
      IFSTSV = IFST;
      ILSTSV = ILST;
      IFST1 = IFST;
      ILST1 = ILST;
      IFST2 = IFST;
      ILST2 = ILST;
      RES = ZERO;

      // Test without accumulating Q

      slaset('Full', N, N, ZERO, ONE, Q, LDT );
      strexc('N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1 );
      for (I = 1; I <= N; I++) { // 40
         for (J = 1; J <= N; J++) { // 30
            if( I == J && Q( I, J ) != ONE ) RES = RES + ONE / EPS;
            IF( I != J && Q( I, J ) != ZERO ) RES = RES + ONE / EPS;
         } // 30
      } // 40

      // Test with accumulating Q

      slaset('Full', N, N, ZERO, ONE, Q, LDT );
      strexc('V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2 );

      // Compare T1 with T2

      for (I = 1; I <= N; I++) { // 60
         for (J = 1; J <= N; J++) { // 50
            if( T1( I, J ) != T2( I, J ) ) RES = RES + ONE / EPS;
         } // 50
      } // 60
      if (IFST1 != IFST2) RES = RES + ONE / EPS;
      if( ILST1 != ILST2 ) RES = RES + ONE / EPS;
      IF( INFO1 != INFO2 ) RES = RES + ONE / EPS;

      // Test for successful reordering of T2

      if ( INFO2 != 0 ) {
         NINFO[INFO2] = NINFO( INFO2 ) + 1;
      } else {
         if( ( IFST2-IFSTSV ).abs() > 1 ) RES = RES + ONE / EPS;
         IF( ( ILST2-ILSTSV ).abs() > 1 ) RES = RES + ONE / EPS;
      }

      // Test for small residual, and orthogonality of Q

      shst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RESULT );
      RES = RES + RESULT( 1 ) + RESULT( 2 );

      // Test for T2 being in Schur form

      LOC = 1;
      } // 70
      if ( T2( LOC+1, LOC ) != ZERO ) {

         // 2 by 2 block

         if( T2( LOC, LOC+1 ) == ZERO || T2( LOC, LOC ) != T2( LOC+1, LOC+1 ) || sign( ONE, T2( LOC, LOC+1 ) ) == sign( ONE, T2( LOC+1, LOC ) ) )RES = RES + ONE / EPS;
         for (I = LOC + 2; I <= N; I++) { // 80
            if( T2( I, LOC ) != ZERO ) RES = RES + ONE / RES;
            IF( T2( I, LOC+1 ) != ZERO ) RES = RES + ONE / RES;
         } // 80
         LOC = LOC + 2;
      } else {

         // 1 by 1 block

         for (I = LOC + 1; I <= N; I++) { // 90
            if( T2( I, LOC ) != ZERO ) RES = RES + ONE / RES;
         } // 90
         LOC = LOC + 1;
      }
      if (LOC < N) GO TO 70;
      if ( RES > RMAX ) {
         RMAX = RES;
         LMAX = KNT;
      }
      GO TO 10;
      }
