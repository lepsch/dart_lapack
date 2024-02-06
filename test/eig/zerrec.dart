import 'package:lapack/src/complex.dart';

import 'common.dart';

      void zerrec(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 4, LW = NMAX*( NMAX+2 ) ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IFST, ILST, INFO, J, M, NT;
      double             SCALE;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      double             RW( LW ), S( NMAX ), SEP( NMAX ), SWORK( NMAX );
      Complex         A( NMAX, NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX ), WORK( LW ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZTREXC, ZTRSEN, ZTRSNA, ZTRSYL
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / infoc / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / srnamc / srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      infoc.OK = true;
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

      // Test ZTRSYL

      srnamc.SRNAMT = 'ZTRSYL';
      infoc.INFOT = 1;
      ztrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      ztrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO );
      chkxer('ZTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      NT = NT + 8;

      // Test ZTRSYL3

      srnamc.SRNAMT = 'ZTRSYL3';
      infoc.INFOT = 1;
      ztrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      ztrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, SWORK, NMAX, INFO );
      chkxer('ZTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      NT = NT + 8;

      // Test ZTREXC

      srnamc.SRNAMT = 'ZTREXC';
      IFST = 1;
      ILST = 1;
      infoc.INFOT = 1;
      ztrexc('X', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrexc('N', -1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ILST = 2;
      ztrexc('N', 2, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztrexc('V', 2, A, 2, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      IFST = 0;
      ILST = 1;
      ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      IFST = 2;
      ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      IFST = 1;
      ILST = 0;
      ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ILST = 2;
      ztrexc('V', 1, A, 1, B, 1, IFST, ILST, INFO );
      chkxer('ZTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      NT = NT + 8;

      // Test ZTRSNA

      srnamc.SRNAMT = 'ZTRSNA';
      infoc.INFOT = 1;
      ztrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      ztrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK, 2, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      ztrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      ztrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 16;
      ztrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK, 1, RW, INFO );
      chkxer('ZTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      NT = NT + 9;

      // Test ZTRSEN

      SEL[1] = false;
      srnamc.SRNAMT = 'ZTRSEN';
      infoc.INFOT = 1;
      ztrsen('X', 'N', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztrsen('N', 'X', SEL, 0, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztrsen('N', 'N', SEL, -1, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztrsen('N', 'N', SEL, 2, A, 1, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 2, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztrsen('N', 'V', SEL, 2, A, 2, B, 1, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 14;
      ztrsen('N', 'V', SEL, 2, A, 2, B, 2, X, M, S( 1 ), SEP( 1 ), WORK, 0, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 14;
      ztrsen('E', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ), WORK, 1, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 14;
      ztrsen('V', 'V', SEL, 3, A, 3, B, 3, X, M, S( 1 ), SEP( 1 ), WORK, 3, INFO );
      chkxer('ZTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      NT = NT + 8;

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( infoc.NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( infoc.NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );
      return;
      }
