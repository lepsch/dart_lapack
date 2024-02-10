      void cerrec(final int PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX, LW;
      const              NMAX = 4, LW = NMAX*( NMAX+2 ) ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, IFST, ILST, INFO, J, M, NT;
      double               SCALE;
      bool               SEL( NMAX );
      double               RW( LW ), S( NMAX ), SEP( NMAX ), SWORK( NMAX );
      Complex            A( NMAX, NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX ), WORK( LW ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, CTREXC, CTRSEN, CTRSNA, CTRSYL, CTRSYL3
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      NOUT = NUNIT;
      OK = true;
      NT = 0;

      // Initialize A, B and SEL

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = ZERO;
            B[I][J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A[I][I] = ONE;
         SEL[I] = true;
      } // 30

      // Test CTRSYL

     srnamc.SRNAMT = 'CTRSYL';
      INFOT = 1;
      ctrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO );
      chkxer('CTRSYL', INFOT, NOUT, LERR, OK );
      NT = NT + 8;

      // Test CTRSYL3

     srnamc.SRNAMT = 'CTRSYL3';
      INFOT = 1;
      ctrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('CTRSYL3', INFOT, NOUT, LERR, OK );
      NT = NT + 8;

      // Test CTREXC

     srnamc.SRNAMT = 'CTREXC';
      IFST = 1;
      ILST = 1;
      INFOT = 1;
      ctrexc('X', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrexc('N', -1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ILST = 2;
      ctrexc('N', 2, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrexc('V', 2, A, 2, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      IFST = 0;
      ILST = 1;
      ctrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      IFST = 2;
      ctrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      IFST = 1;
      ILST = 0;
      ctrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ILST = 2;
      ctrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('CTREXC', INFOT, NOUT, LERR, OK );
      NT = NT + 8;

      // Test CTRSNA

     srnamc.SRNAMT = 'CTRSNA';
      INFOT = 1;
      ctrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ctrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ctrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ctrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ctrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 16;
      ctrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK, 1, RW, INFO );
      chkxer('CTRSNA', INFOT, NOUT, LERR, OK );
      NT = NT + 9;

      // Test CTRSEN

      SEL[1] = false;
     srnamc.SRNAMT = 'CTRSEN';
      INFOT = 1;
      ctrsen('X', 'N', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrsen('N', 'X', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrsen('N', 'N', SEL, -1, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsen('N', 'N', SEL, 2, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 2, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ctrsen('N', 'V', SEL, 2, A, 2, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 14;
      ctrsen('N', 'V', SEL, 2, A, 2, B, 2, X, M, S( 1 ), SEP( 1 ), WORK, 0, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 14;
      ctrsen('E', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 14;
      ctrsen('V', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ), WORK, 3, INFO );
      chkxer('CTRSEN', INFOT, NOUT, LERR, OK );
      NT = NT + 8;

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT(' ${.a3} routines passed the tests of the error exits (${.i3} tests done)' );
 9998 FORMAT( ' *** ${.a3} routines failed the tests of the error exits ***' );
      }
