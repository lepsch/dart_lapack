      SUBROUTINE SCHKEC( THRESH, TSTERR, NIN, NOUT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NIN, NOUT;
      REAL               THRESH;
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               OK;
      String             PATH;
      int                KLAEXC, KLALN2, KLANV2, KLAQTR, KLASY2, KTREXC, KTRSEN, KTRSNA, KTRSYL, KTRSYL3, LLAEXC, LLALN2, LLANV2, LLAQTR, LLASY2, LTREXC, LTRSYL, NLANV2, NLAQTR, NLASY2, NTESTS, NTRSYL, KTGEXC, LTGEXC;
      REAL               EPS, RLAEXC, RLALN2, RLANV2, RLAQTR, RLASY2, RTREXC, SFMIN, RTGEXC;
      // ..
      // .. Local Arrays ..
      int                FTRSYL( 3 ), ITRSYL( 2 ), LTRSEN( 3 ), LTRSNA( 3 ), NLAEXC( 2 ), NLALN2( 2 ), NTGEXC( 2 ), NTREXC( 3 ), NTRSEN( 3 ), NTRSNA( 3 );
      REAL               RTRSEN( 3 ), RTRSNA( 3 ), RTRSYL( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SERREC, SGET31, SGET32, SGET33, SGET34, SGET35, SGET36, SGET37, SGET38, SGET39, SGET40, SSYL01
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Single precision';
      PATH( 2: 3 ) = 'EC';
      EPS = SLAMCH( 'P' );
      SFMIN = SLAMCH( 'S' );

      // Print header information

      WRITE( NOUT, FMT = 9989 );
      WRITE( NOUT, FMT = 9988 )EPS, SFMIN;
      WRITE( NOUT, FMT = 9987 )THRESH;

      // Test error exits if TSTERR is true;

      if (TSTERR) CALL SERREC( PATH, NOUT );

      OK = true;
      sget31(RLALN2, LLALN2, NLALN2, KLALN2 );
      if ( RLALN2 > THRESH || NLALN2( 1 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9999 )RLALN2, LLALN2, NLALN2, KLALN2;
      }

      sget32(RLASY2, LLASY2, NLASY2, KLASY2 );
      if ( RLASY2 > THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9998 )RLASY2, LLASY2, NLASY2, KLASY2;
      }

      sget33(RLANV2, LLANV2, NLANV2, KLANV2 );
      if ( RLANV2 > THRESH || NLANV2 != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9997 )RLANV2, LLANV2, NLANV2, KLANV2;
      }

      sget34(RLAEXC, LLAEXC, NLAEXC, KLAEXC );
      if ( RLAEXC > THRESH || NLAEXC( 2 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9996 )RLAEXC, LLAEXC, NLAEXC, KLAEXC;
      }

      sget35(RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL );
      if ( RTRSYL( 1 ) > THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9995 )RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL;
      }

      ssyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3 );
      if ( FTRSYL( 1 ) > 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9970 )FTRSYL( 1 ), RTRSYL( 1 ), THRESH;
      }
      if ( FTRSYL( 2 ) > 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9971 )FTRSYL( 2 ), RTRSYL( 2 ), THRESH;
      }
      if ( FTRSYL( 3 ) > 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9972 )FTRSYL( 3 );
      }

      sget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN );
      if ( RTREXC > THRESH || NTREXC( 3 ) > 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9994 )RTREXC, LTREXC, NTREXC, KTREXC;
      }

      sget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN );
      if ( RTRSNA( 1 ) > THRESH || RTRSNA( 2 ) > THRESH || NTRSNA( 1 ) != 0 || NTRSNA( 2 ) != 0 || NTRSNA( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9993 )RTRSNA, LTRSNA, NTRSNA, KTRSNA;
      }

      sget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN );
      if ( RTRSEN( 1 ) > THRESH || RTRSEN( 2 ) > THRESH || NTRSEN( 1 ) != 0 || NTRSEN( 2 ) != 0 || NTRSEN( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9992 )RTRSEN, LTRSEN, NTRSEN, KTRSEN;
      }

      sget39(RLAQTR, LLAQTR, NLAQTR, KLAQTR );
      if ( RLAQTR > THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9991 )RLAQTR, LLAQTR, NLAQTR, KLAQTR;
      }

      sget40(RTGEXC, LTGEXC, NTGEXC, KTGEXC, NIN );
      if ( RTGEXC > THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9986 )RTGEXC, LTGEXC, NTGEXC, KTGEXC;
      }

      NTESTS = KLALN2 + KLASY2 + KLANV2 + KLAEXC + KTRSYL + KTREXC + KTRSNA + KTRSEN + KLAQTR       IF( OK ) WRITE( NOUT, FMT = 9990 )PATH, NTESTS;

      return;
 9999 FORMAT( ' Error in SLALN2: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 );
 9998 FORMAT( ' Error in SLASY2: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 );
 9997 FORMAT( ' Error in SLANV2: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 );
 9996 FORMAT( ' Error in SLAEXC: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 );
 9995 FORMAT( ' Error in STRSYL: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 );
 9994 FORMAT( ' Error in STREXC: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', 3I8, ' KNT=', I8 );
 9993 FORMAT( ' Error in STRSNA: RMAX =', 3E12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 );
 9992 FORMAT( ' Error in STRSEN: RMAX =', 3E12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=', I8 );
 9991 FORMAT( ' Error in SLAQTR: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', I8, ' KNT=', I8 );
 9990 FORMAT( / 1X, 'All tests for ', A3, ' routines passed the thresh', 'old ( ', I6, ' tests run)' );
 9989 FORMAT( ' Tests of the Nonsymmetric eigenproblem condition estim', 'ation routines', / ' SLALN2, SLASY2, SLANV2, SLAEXC, STRS', 'YL, STREXC, STRSNA, STRSEN, SLAQTR', / );
 9988 FORMAT( ' Relative machine precision (EPS) = ', E16.6, / ' Safe ', 'minimum (SFMIN)             = ', E16.6, / );
 9987 FORMAT( ' Routines pass computational tests if test ratio is les', 's than', F8.2, / / );
 9986 FORMAT( ' Error in STGEXC: RMAX =', E12.3, / ' LMAX = ', I8, ' N', 'INFO=', 2I8, ' KNT=', I8 );
 9972 FORMAT( 'STRSYL and STRSYL3 compute an inconsistent result ', 'factor in ', I8, ' tests.');
 9971 FORMAT( 'Error in STRSYL3: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 );
 9970 FORMAT( 'Error in STRSYL: ', I8, ' tests fail the threshold.', / 'Maximum test ratio =', D12.3, ' threshold =', D12.3 );

      // End of SCHKEC

      }
