      SUBROUTINE SERREC( PATH, NUNIT )

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
      REAL               ONE, ZERO
      const              NMAX = 4, ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, IFST, ILST, INFO, J, M, NT;
      REAL               SCALE
      // ..
      // .. Local Arrays ..
      bool               SEL( NMAX );
      int                IWORK( NMAX );
      REAL               A( NMAX, NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX ), S( NMAX ), SEP( NMAX ), WI( NMAX ), WORK( NMAX ), WR( NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, STREXC, STRSEN, STRSNA, STRSYL, STRSYL3
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
   10    CONTINUE
   20 CONTINUE
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE
         SEL( I ) = .TRUE.
   30 CONTINUE

      // Test STRSYL

      SRNAMT = 'STRSYL'
      INFOT = 1
      strsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 2
      strsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 3
      strsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 4
      strsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      strsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 7
      strsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 9
      strsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      INFOT = 11
      strsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO );
      chkxer('STRSYL', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test STRSYL3

      SRNAMT = 'STRSYL3'
      INFOT = 1
      strsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 2
      strsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 3
      strsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 4
      strsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 5
      strsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 7
      strsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 9
      strsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      INFOT = 11
      strsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX, INFO );
      chkxer('STRSYL3', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test STREXC

      SRNAMT = 'STREXC'
      IFST = 1
      ILST = 1
      INFOT = 1
      strexc('X', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 2
      strexc('N', -1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ILST = 2
      strexc('N', 2, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 6
      strexc('V', 2, A, 2, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7
      IFST = 0
      ILST = 1
      strexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 7
      IFST = 2
      strexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8
      IFST = 1
      ILST = 0
      strexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      INFOT = 8
      ILST = 2
      strexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO );
      chkxer('STREXC', INFOT, NOUT, LERR, OK );
      NT = NT + 8

      // Test STRSNA

      SRNAMT = 'STRSNA'
      INFOT = 1
      strsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 2
      strsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 4
      strsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK, 1, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 6
      strsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 8
      strsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 10
      strsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK, 2, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13
      strsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK, 1, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 13
      strsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK, 2, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      INFOT = 16
      strsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK, 1, IWORK, INFO );
      chkxer('STRSNA', INFOT, NOUT, LERR, OK );
      NT = NT + 9

      // Test STRSEN

      SEL( 1 ) = .FALSE.
      SRNAMT = 'STRSEN'
      INFOT = 1
      strsen('X', 'N', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 2
      strsen('N', 'X', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 4
      strsen('N', 'N', SEL, -1, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 6
      strsen('N', 'N', SEL, 2, A, 1, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 2, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 8
      strsen('N', 'V', SEL, 2, A, 2, B, 1, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      strsen('N', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 0, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      strsen('E', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 15
      strsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 3, IWORK, 2, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 17
      strsen('E', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 1, IWORK, 0, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
      INFOT = 17
      strsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ), SEP( 1 ), WORK, 4, IWORK, 1, INFO );
      chkxer('STRSEN', INFOT, NOUT, LERR, OK );
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

      // End of SERREC

      }
