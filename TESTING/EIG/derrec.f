      SUBROUTINE DERREC( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      double             ONE, ZERO;
      const              NMAX = 4, ONE = 1.0D0, ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, IFST, ILST, INFO, J, M, NT;
      double             SCALE;
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IWORK( NMAX );
      double             A( NMAX, NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX ), S( NMAX ), SEP( NMAX ), WI( NMAX ), WORK( NMAX ), WR( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DTREXC, DTRSEN, DTRSNA, DTRSYL, DTRSYL3
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      OK = .TRUE.
      NT = 0

      // Initialize A, B and SEL

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = ZERO
            B( I, J ) = ZERO
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE
         SEL( I ) = .TRUE.
      } // 30

      // Test DTRSYL

      SRNAMT = 'DTRSYL'
      INFOT = 1
      dtrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dtrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dtrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO );
      chkxer('DTRSYL', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test DTRSYL3

      SRNAMT = 'DTRSYL3'
      INFOT = 1
      dtrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dtrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dtrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('DTRSYL3', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test DTREXC

      SRNAMT = 'DTREXC'
      IFST = 1
      ILST = 1
      INFOT = 1
      dtrexc('X', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrexc('N', -1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ILST = 2
      dtrexc('N', 2, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrexc('V', 2, A, 2, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7
      IFST = 0
      ILST = 1
      dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7
      IFST = 2
      dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8
      IFST = 1
      ILST = 0
      dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8
      ILST = 2
      dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('DTREXC', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test DTRSNA

      SRNAMT = 'DTRSNA'
      INFOT = 1
      dtrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dtrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13
      dtrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK, 1, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13
      dtrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK, 2, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 16
      dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK, 1, IWORK, INFO );
      chkxer('DTRSNA', INFOT, NOUT, LERR, OK );
      NT = NT + 9

      // Test DTRSEN

      SEL( 1 ) = .FALSE.
      SRNAMT = 'DTRSEN'
      INFOT = 1
      dtrsen('X', 'N', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrsen('N', 'X', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrsen('N', 'N', SEL, -1, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsen('N', 'N', SEL, 2, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 2, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dtrsen('N', 'V', SEL, 2, A, 2, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      dtrsen('N', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 0, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      dtrsen('E', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 3, IWORK, 2, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 17
      dtrsen('E', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 0, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 17
      dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 4, IWORK, 1, INFO );
      chkxer('DTRSEN', INFOT, NOUT, LERR, OK );
      NT = NT + 10

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

      RETURN
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ex', 'its ***' )

      // End of DERREC

      }
