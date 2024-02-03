      SUBROUTINE CCHKEC( THRESH, TSTERR, NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NIN, NOUT;
      REAL               THRESH
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               OK;
      String             PATH;
      int                KTREXC, KTRSEN, KTRSNA, KTRSYL, KTRSYL3, LTREXC, LTRSYL, NTESTS, NTREXC, NTRSYL;
      REAL               EPS, RTREXC, SFMIN
      // ..
      // .. Local Arrays ..
      int                FTRSYL( 3 ), ITRSYL( 2 ), LTRSEN( 3 ), LTRSNA( 3 ), NTRSEN( 3 ), NTRSNA( 3 );
      REAL               RTRSEN( 3 ), RTRSNA( 3 ), RTRSYL( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CERREC, CGET35, CGET36, CGET37, CGET38, CSYL01
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'EC'
      EPS = SLAMCH( 'P' )
      SFMIN = SLAMCH( 'S' )
      WRITE( NOUT, FMT = 9994 )
      WRITE( NOUT, FMT = 9993 )EPS, SFMIN
      WRITE( NOUT, FMT = 9992 )THRESH

      // Test error exits if TSTERR is .TRUE.

      IF( TSTERR ) CALL CERREC( PATH, NOUT )

      OK = .TRUE.
      cget35(RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL, NIN );
      if ( RTRSYL( 1 ).GT.THRESH ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9999 )RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL
      }

      csyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3 );
      if ( FTRSYL( 1 ).GT.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9970 )FTRSYL( 1 ), RTRSYL( 1 ), THRESH
      }
      if ( FTRSYL( 2 ).GT.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9971 )FTRSYL( 2 ), RTRSYL( 2 ), THRESH
      }
      if ( FTRSYL( 3 ).GT.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9972 )FTRSYL( 3 )
      }

      cget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN );
      if ( RTREXC.GT.THRESH .OR. NTREXC.GT.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9998 )RTREXC, LTREXC, NTREXC, KTREXC
      }

      cget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN );
      if ( RTRSNA( 1 ).GT.THRESH .OR. RTRSNA( 2 ).GT.THRESH .OR. NTRSNA( 1 ).NE.0 .OR. NTRSNA( 2 ).NE.0 .OR. NTRSNA( 3 ).NE.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9997 )RTRSNA, LTRSNA, NTRSNA, KTRSNA
      }

      cget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN );
      if ( RTRSEN( 1 ).GT.THRESH .OR. RTRSEN( 2 ).GT.THRESH .OR. NTRSEN( 1 ).NE.0 .OR. NTRSEN( 2 ).NE.0 .OR. NTRSEN( 3 ).NE.0 ) {
         OK = .FALSE.
         WRITE( NOUT, FMT = 9996 )RTRSEN, LTRSEN, NTRSEN, KTRSEN
      }

      NTESTS = KTRSYL + KTREXC + KTRSNA + KTRSEN
      IF( OK ) WRITE( NOUT, FMT = 9995 )PATH, NTESTS

 9999 FORMAT( ' Error in CTRSYL: RMAX =', E12.3, / ' LMAX = ', I8, ' NINFO=', I8, ' KNT=', I8 )
 9998 FORMAT( ' Error in CTREXC: RMAX =', E12.3, / ' LMAX = ', I8, ' NINFO=', I8, ' KNT=', I8 )
 9997 FORMAT( ' Error in CTRSNA: RMAX =', 3E12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 )
 9996 FORMAT( ' Error in CTRSEN: RMAX =', 3E12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 )
 9995 FORMAT( / 1X, 'All tests for ', A3, ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( ' Tests of the Nonsymmetric eigenproblem condition', ' estimation routines', / ' CTRSYL, CTREXC, CTRSNA, CTRSEN', / )
 9993 FORMAT( ' Relative machine precision (EPS) = ', E16.6, / ' Safe minimum (SFMIN)             = ', E16.6, / )
 9992 FORMAT( ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / / )
 9972 FORMAT( 'CTRSYL and CTRSYL3 compute an inconsistent scale ', 'factor in ', I8, ' tests.')
 9971 FORMAT( 'Error in CTRSYL3: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 )
 9970 FORMAT( 'Error in CTRSYL: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 )
      RETURN

      // End of CCHKEC

      }
