      void cget36(final int RMAX, final int LMAX, final int NINFO, final int KNT, final int NIN,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KNT, LMAX, NIN, NINFO;
      double               RMAX;
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                LDT, LWORK;
      const              LDT = 10, LWORK = 2*LDT*LDT ;
      int                I, IFST, ILST, INFO1, INFO2, J, N;
      double               EPS, RES;
      Complex            CTEMP;
      double               RESULT( 2 ), RWORK( LDT );
      Complex            DIAG( LDT ), Q( LDT, LDT ), T1( LDT, LDT ), T2( LDT, LDT ), TMP( LDT, LDT ), WORK( LWORK );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHST01, CLACPY, CLASET, CTREXC

      EPS = SLAMCH( 'P' );
      RMAX = ZERO;
      LMAX = 0;
      KNT = 0;
      NINFO = 0;

      // Read input data until N=0

      } // 10
      READ( NIN, FMT = * )N, IFST, ILST;
      if (N == 0) return;
      KNT = KNT + 1;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N );
      } // 20
      clacpy('F', N, N, TMP, LDT, T1, LDT );
      clacpy('F', N, N, TMP, LDT, T2, LDT );
      RES = ZERO;

      // Test without accumulating Q

      claset('Full', N, N, CZERO, CONE, Q, LDT );
      ctrexc('N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1 );
      for (I = 1; I <= N; I++) { // 40
         for (J = 1; J <= N; J++) { // 30
            if( I == J && Q( I, J ) != CONE ) RES = RES + ONE / EPS;
            IF( I != J && Q( I, J ) != CZERO ) RES = RES + ONE / EPS;
         } // 30
      } // 40

      // Test with accumulating Q

      claset('Full', N, N, CZERO, CONE, Q, LDT );
      ctrexc('V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2 );

      // Compare T1 with T2

      for (I = 1; I <= N; I++) { // 60
         for (J = 1; J <= N; J++) { // 50
            if( T1( I, J ) != T2( I, J ) ) RES = RES + ONE / EPS;
         } // 50
      } // 60
      if (INFO1 != 0 || INFO2 != 0) NINFO = NINFO + 1;
      IF( INFO1 != INFO2 ) RES = RES + ONE / EPS;

      // Test for successful reordering of T2

      ccopy(N, TMP, LDT+1, DIAG, 1 );
      if ( IFST < ILST ) {
         for (I = IFST + 1; I <= ILST; I++) { // 70
            CTEMP = DIAG( I );
            DIAG[I] = DIAG( I-1 );
            DIAG[I-1] = CTEMP;
         } // 70
      } else if ( IFST > ILST ) {
         for (I = IFST - 1; I >= ILST; I--) { // 80
            CTEMP = DIAG( I+1 );
            DIAG[I+1] = DIAG( I );
            DIAG[I] = CTEMP;
         } // 80
      }
      for (I = 1; I <= N; I++) { // 90
         if( T2( I, I ) != DIAG( I ) ) RES = RES + ONE / EPS;
      } // 90

      // Test for small residual, and orthogonality of Q

      chst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT );
      RES = RES + RESULT( 1 ) + RESULT( 2 );

      // Test for T2 being in Schur form

      for (J = 1; J <= N - 1; J++) { // 110
         for (I = J + 1; I <= N; I++) { // 100
            if( T2( I, J ) != CZERO ) RES = RES + ONE / EPS;
         } // 100
      } // 110
      if ( RES > RMAX ) {
         RMAX = RES;
         LMAX = KNT;
      }
      GO TO 10;
      }