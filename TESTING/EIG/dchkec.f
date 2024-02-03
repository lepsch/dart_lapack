      SUBROUTINE DCHKEC( THRESH, TSTERR, NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NIN, NOUT;
      double             THRESH;
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               OK;
      String             PATH;
      int                KLAEXC, KLALN2, KLANV2, KLAQTR, KLASY2, KTREXC, KTRSEN, KTRSNA, KTRSYL, KTRSYL3, LLAEXC, LLALN2, LLANV2, LLAQTR, LLASY2, LTREXC, LTRSYL, NLANV2, NLAQTR, NLASY2, NTESTS, NTRSYL, KTGEXC, LTGEXC;
      double             EPS, RLAEXC, RLALN2, RLANV2, RLAQTR, RLASY2, RTREXC, SFMIN, RTGEXC;
      // ..
      // .. Local Arrays ..
      int                FTRSYL( 3 ), ITRSYL( 2 ), LTRSEN( 3 ), LTRSNA( 3 ), NLAEXC( 2 ), NLALN2( 2 ), NTGEXC( 2 ), NTREXC( 3 ), NTRSEN( 3 ), NTRSNA( 3 );
      double             RTRSEN( 3 ), RTRSNA( 3 ), RTRSYL( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DERREC, DGET31, DGET32, DGET33, DGET34, DGET35, DGET36, DGET37, DGET38, DGET39, DGET40, DSYL01
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'EC'
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )

      // Print header information

      WRITE( NOUT, FMT = 9989 )
      WRITE( NOUT, FMT = 9988 )EPS, SFMIN
      WRITE( NOUT, FMT = 9987 )THRESH

      // Test error exits if TSTERR is true;

      if (TSTERR) CALL DERREC( PATH, NOUT );

      OK = true;
      dget31(RLALN2, LLALN2, NLALN2, KLALN2 );
      if ( RLALN2.GT.THRESH || NLALN2( 1 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9999 )RLALN2, LLALN2, NLALN2, KLALN2
      }

      dget32(RLASY2, LLASY2, NLASY2, KLASY2 );
      if ( RLASY2.GT.THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9998 )RLASY2, LLASY2, NLASY2, KLASY2
      }

      dget33(RLANV2, LLANV2, NLANV2, KLANV2 );
      if ( RLANV2.GT.THRESH || NLANV2 != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9997 )RLANV2, LLANV2, NLANV2, KLANV2
      }

      dget34(RLAEXC, LLAEXC, NLAEXC, KLAEXC );
      if ( RLAEXC.GT.THRESH || NLAEXC( 2 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9996 )RLAEXC, LLAEXC, NLAEXC, KLAEXC
      }

      dget35(RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL );
      if ( RTRSYL( 1 ).GT.THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9995 )RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL
      }

      dsyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3 );
      if ( FTRSYL( 1 ).GT.0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9970 )FTRSYL( 1 ), RTRSYL( 1 ), THRESH
      }
      if ( FTRSYL( 2 ).GT.0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9971 )FTRSYL( 2 ), RTRSYL( 2 ), THRESH
      }
      if ( FTRSYL( 3 ).GT.0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9972 )FTRSYL( 3 )
      }

      dget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN );
      if ( RTREXC.GT.THRESH || NTREXC( 3 ).GT.0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9994 )RTREXC, LTREXC, NTREXC, KTREXC
      }

      dget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN );
      if ( RTRSNA( 1 ).GT.THRESH || RTRSNA( 2 ).GT.THRESH || NTRSNA( 1 ) != 0 || NTRSNA( 2 ) != 0 || NTRSNA( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9993 )RTRSNA, LTRSNA, NTRSNA, KTRSNA
      }

      dget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN );
      if ( RTRSEN( 1 ).GT.THRESH || RTRSEN( 2 ).GT.THRESH || NTRSEN( 1 ) != 0 || NTRSEN( 2 ) != 0 || NTRSEN( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9992 )RTRSEN, LTRSEN, NTRSEN, KTRSEN
      }

      dget39(RLAQTR, LLAQTR, NLAQTR, KLAQTR );
      if ( RLAQTR.GT.THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9991 )RLAQTR, LLAQTR, NLAQTR, KLAQTR
      }

      dget40(RTGEXC, LTGEXC, NTGEXC, KTGEXC, NIN );
      if ( RTGEXC.GT.THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9986 )RTGEXC, LTGEXC, NTGEXC, KTGEXC
      }

      NTESTS = KLALN2 + KLASY2 + KLANV2 + KLAEXC + KTRSYL + KTREXC + KTRSNA + KTRSEN + KLAQTR + KTGEXC       IF( OK ) WRITE( NOUT, FMT = 9990 )PATH, NTESTS

      RETURN
 9999 FORMAT( ' Error in DLALN2: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 )
 9998 FORMAT( ' Error in DLASY2: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 )
 9997 FORMAT( ' Error in DLANV2: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 )
 9996 FORMAT( ' Error in DLAEXC: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 )
 9995 FORMAT( ' Error in DTRSYL: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 )
 9994 FORMAT( ' Error in DTREXC: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', 3I8, ' KNT=', I8 )
 9993 FORMAT( ' Error in DTRSNA: RMAX =', 3D12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 )
 9992 FORMAT( ' Error in DTRSEN: RMAX =', 3D12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 )
 9991 FORMAT( ' Error in DLAQTR: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 )
 9990 FORMAT( / 1X, 'All tests for ', A3, ' routines passed the thresh', 'old ( ', I6, ' tests run)' )
 9989 FORMAT( ' Tests of the Nonsymmetric eigenproblem condition estim', 'ation routines', / ' DLALN2, DLASY2, DLANV2, DLAEXC, DTRS', 'YL, DTREXC, DTRSNA, DTRSEN, DLAQTR, DTGEXC', / )
 9988 FORMAT( ' Relative machine precision (EPS) = ', D16.6, / ' Safe ', 'minimum (SFMIN)             = ', D16.6, / )
 9987 FORMAT( ' Routines pass computational tests if test ratio is les', 's than', F8.2, / / )
 9986 FORMAT( ' Error in DTGEXC: RMAX =', D12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 )
 9972 FORMAT( 'DTRSYL and DTRSYL3 compute an inconsistent result ', 'factor in ', I8, ' tests.')
 9971 FORMAT( 'Error in DTRSYL3: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 )
 9970 FORMAT( 'Error in DTRSYL: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 )

      // End of DCHKEC

      }
