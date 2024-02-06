      void main() {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                MAXIN;
      const              MAXIN = 12 ;
      int                NMAX;
      const              NMAX =  50 ;
      int                MAXRHS;
      const              MAXRHS = 16 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      bool               FATAL, TSTERR;
      int                VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      int                I, NN, NNS, NNT;
      double             EPS, S1, S2, THRESH;

      int                NVAL( MAXIN ), NSVAL( MAXIN ), NTVAL( NTYPES );
      Complex         WORKA( NMAX, NMAX );
      Complex         WORKASAV( NMAX, NMAX );
      Complex         WORKB( NMAX, MAXRHS );
      Complex         WORKXACT( NMAX, MAXRHS );
      Complex         WORKBSAV( NMAX, MAXRHS );
      Complex         WORKX( NMAX, MAXRHS );
      Complex         WORKAFAC( NMAX, NMAX );
      Complex         WORKAINV( NMAX, NMAX );
      Complex         WORKARF( (NMAX*(NMAX+1))/2 );
      Complex         WORKAP( (NMAX*(NMAX+1))/2 );
      Complex         WORKARFINV( (NMAX*(NMAX+1))/2 );
      Complex         Z_WORK_ZLATMS( 3 * NMAX );
      Complex         Z_WORK_ZPOT02( NMAX, MAXRHS );
      Complex         Z_WORK_ZPOT03( NMAX, NMAX );
      double             D_WORK_ZLATMS( NMAX );
      double             D_WORK_ZLANHE( NMAX );
      double             D_WORK_ZPOT01( NMAX );
      double             D_WORK_ZPOT02( NMAX );
      double             D_WORK_ZPOT03( NMAX );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DSECND;
      // EXTERNAL DLAMCH, DSECND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ILAVER, ZDRVRFP, ZDRVRF1, ZDRVRF2, ZDRVRF3, ZDRVRF4

      S1 = DSECND( );
      FATAL = false;

      // Read a dummy line.

      READ( NIN, FMT = * );

      // Report LAPACK version tag (e.g. LAPACK-3.2.0)

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH;

      // Read the values of N

      READ( NIN, FMT = * )NN;
      if ( NN < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1;
         NN = 0;
         FATAL = true;
      } else if ( NN > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN;
         NN = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN );
      for (I = 1; I <= NN; I++) { // 10
         if ( NVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )' M  ', NVAL( I ), 0;
            FATAL = true;
         } else if ( NVAL( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9995 )' M  ', NVAL( I ), NMAX;
            FATAL = true;
         }
      } // 10
      if (NN > 0) WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN );

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS;
      if ( NNS < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1;
         NNS = 0;
         FATAL = true;
      } else if ( NNS > MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN;
         NNS = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS );
      for (I = 1; I <= NNS; I++) { // 30
         if ( NSVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0;
            FATAL = true;
         } else if ( NSVAL( I ) > MAXRHS ) {
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS;
            FATAL = true;
         }
      } // 30
      if (NNS > 0) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS );

      // Read the matrix types

      READ( NIN, FMT = * )NNT;
      if ( NNT < 1 ) {
         WRITE( NOUT, FMT = 9996 )' NMA', NNT, 1;
         NNT = 0;
         FATAL = true;
      } else if ( NNT > NTYPES ) {
         WRITE( NOUT, FMT = 9995 )' NMA', NNT, NTYPES;
         NNT = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NTVAL( I ), I = 1, NNT );
      for (I = 1; I <= NNT; I++) { // 320
         if ( NTVAL( I ) < 0 ) {
            WRITE( NOUT, FMT = 9996 )'TYPE', NTVAL( I ), 0;
            FATAL = true;
         } else if ( NTVAL( I ) > NTYPES ) {
            WRITE( NOUT, FMT = 9995 )'TYPE', NTVAL( I ), NTYPES;
            FATAL = true;
         }
      } // 320
      if (NNT > 0) WRITE( NOUT, FMT = 9993 )'TYPE', ( NTVAL( I ), I = 1, NNT );

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH;
      WRITE( NOUT, FMT = 9992 )THRESH;

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR;

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 );
         STOP;
      }

      // Calculate and print the machine dependent constants.

      EPS = dlamch( 'Underflow threshold' );
      WRITE( NOUT, FMT = 9991 )'underflow', EPS;
      EPS = dlamch( 'Overflow threshold' );
      WRITE( NOUT, FMT = 9991 )'overflow ', EPS;
      EPS = dlamch( 'Epsilon' );
      WRITE( NOUT, FMT = 9991 )'precision', EPS;
      WRITE( NOUT, FMT = * );

      // Test the error exit of:

      if (TSTERR) zerrrfp( NOUT );

// Test the routines: zpftrf, zpftri, zpftrs (as in ZDRVPO).
// This also tests the routines: ztfsm, ztftri, ztfttr, ztrttf.

      zdrvrfp(NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, WORKA, WORKASAV, WORKAFAC, WORKAINV, WORKB, WORKBSAV, WORKXACT, WORKX, WORKARF, WORKARFINV, Z_WORK_ZLATMS, Z_WORK_ZPOT02, Z_WORK_ZPOT03, D_WORK_ZLATMS, D_WORK_ZLANHE, D_WORK_ZPOT01, D_WORK_ZPOT02, D_WORK_ZPOT03 );

// Test the routine: zlanhf

      zdrvrf1(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, D_WORK_ZLANHE );

// Test the conversion routines:
        // zhfttp, ztpthf, ztfttr, ztrttf, ztrttp and ztpttr.

      zdrvrf2(NOUT, NN, NVAL, WORKA, NMAX, WORKARF, WORKAP, WORKASAV );

// Test the routine: ztfsm

      zdrvrf3(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, WORKAINV, WORKAFAC, D_WORK_ZLANHE, Z_WORK_ZPOT03, Z_WORK_ZPOT02 );


// Test the routine: zhfrk

      zdrvrf4(NOUT, NN, NVAL, THRESH, WORKA, WORKAFAC, NMAX, WORKARF, WORKAINV, NMAX,D_WORK_ZLANHE);

      CLOSE ( NIN );
      S2 = DSECND( );
      WRITE( NOUT, FMT = 9998 );
      WRITE( NOUT, FMT = 9997 )S2 - S1;

 9999 FORMAT('\n Execution not attempted due to input errors' );
 9998 FORMAT('\n End of tests' );
 9997 FORMAT( ' Total time used = ${.f12_2} seconds', / );
 9996 FORMAT( ' !! Invalid input value: ${.a4}=${.i6}; must be >=', I6 )
 9995 FORMAT( ' !! Invalid input value: ${.a4}=${.i6}; must be <=', I6 )
 9994 FORMAT('\n Tests of the Complex LAPACK RFP routines \n LAPACK VERSION ${.i1}.${.i1}.', I1, / / ' The following parameter values will be used:' );
 9993 FORMAT('    ${.a4}:  ', 10I6, / 11X, 10I6 );
 9992 FORMAT('\n Routines pass computational tests if test ratio is less than', F8.2, / );
 9991 FORMAT( ' Relative machine ${} is taken to be', D16.6 );
      }
