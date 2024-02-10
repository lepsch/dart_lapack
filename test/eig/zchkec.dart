      void zchkec(final int THRESH, final int TSTERR, final int NIN, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NIN, NOUT;
      double             THRESH;
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               OK;
      String             PATH;
      int                KTREXC, KTRSEN, KTRSNA, KTRSYL, KTRSYL3, LTREXC, LTRSYL, NTESTS, NTREXC, NTRSYL;
      double             EPS, RTREXC, SFMIN;
      int                FTRSYL( 3 ), ITRSYL( 2 ), LTRSEN( 3 ), LTRSNA( 3 ), NTRSEN( 3 ), NTRSNA( 3 );
      double             RTRSEN( 3 ), RTRSNA( 3 ), RTRSYL( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZERREC, ZGET35, ZGET36, ZGET37, ZGET38, ZSYL01
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'EC';
      EPS = dlamch( 'P' );
      SFMIN = dlamch( 'S' );
      WRITE( NOUT, FMT = 9994 );
      WRITE( NOUT, FMT = 9993 )EPS, SFMIN;
      WRITE( NOUT, FMT = 9992 )THRESH;

      // Test error exits if TSTERR is true;

      if (TSTERR) zerrec( PATH, NOUT );

      OK = true;
      zget35(RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL, NIN );
      if ( RTRSYL( 1 ) > THRESH ) {
         OK = false;
         WRITE( NOUT, FMT = 9999 )RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL;
      }

      zsyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3 );
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

      zget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN );
      if ( RTREXC > THRESH || NTREXC > 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9998 )RTREXC, LTREXC, NTREXC, KTREXC;
      }

      zget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN );
      if ( RTRSNA( 1 ) > THRESH || RTRSNA( 2 ) > THRESH || NTRSNA( 1 ) != 0 || NTRSNA( 2 ) != 0 || NTRSNA( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9997 )RTRSNA, LTRSNA, NTRSNA, KTRSNA;
      }

      zget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN );
      if ( RTRSEN( 1 ) > THRESH || RTRSEN( 2 ) > THRESH || NTRSEN( 1 ) != 0 || NTRSEN( 2 ) != 0 || NTRSEN( 3 ) != 0 ) {
         OK = false;
         WRITE( NOUT, FMT = 9996 )RTRSEN, LTRSEN, NTRSEN, KTRSEN;
      }

      NTESTS = KTRSYL + KTRSYL3 + KTREXC + KTRSNA + KTRSEN;
      if (OK) WRITE( NOUT, FMT = 9995 )PATH, NTESTS;

 9999 FORMAT( ' Error in ZTRSYL: RMAX =', D12.3, / ' LMAX = ${.i8} NINFO=${.i8} KNT=${.i8}');
 9998 FORMAT( ' Error in ZTREXC: RMAX =', D12.3, / ' LMAX = ${.i8} NINFO=${.i8} KNT=${.i8}');
 9997 FORMAT( ' Error in ZTRSNA: RMAX =', 3D12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=${.i8}');
 9996 FORMAT( ' Error in ZTRSEN: RMAX =', 3D12.3, / ' LMAX = ', 3I8, ' NINFO=', 3I8, ' KNT=${.i8}');
 9995 FORMAT('\n All tests for ${.a3} routines passed the threshold ( ${.i6} tests run)' );
 9994 FORMAT( ' Tests of the Nonsymmetric eigenproblem condition estimation routines\n ZTRSYL, ZTREXC, ZTRSNA, ZTRSEN\n');
 9993 FORMAT( ' Relative machine precision (EPS) = ', D16.6, / ' Safe minimum (SFMIN)             = ', D16.6, / );
 9992 FORMAT( ' Routines pass computational tests if test ratio is less than', F8.2, / / );
 9970 FORMAT( 'Error in ZTRSYL: ${.i8} tests fail the threshold.\nMaximum test ratio =${.d12_3} threshold =', D12.3 );
 9971 FORMAT( 'Error in ZTRSYL3: ${.i8} tests fail the threshold.\nMaximum test ratio =${.d12_3} threshold =', D12.3 );
 9972 FORMAT( 'ZTRSYL and ZTRSYL3 compute an inconsistent scale factor in ${.i8} tests.');
      }
