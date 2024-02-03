      PROGRAM ZCHKRFP

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
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
      // ..
      // .. Local Scalars ..
      bool               FATAL, TSTERR;
      int                VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      int                I, NN, NNS, NNT;
      double             EPS, S1, S2, THRESH;

      // ..
      // .. Local Arrays ..
      int                NVAL( MAXIN ), NSVAL( MAXIN ), NTVAL( NTYPES );
      COMPLEX*16         WORKA( NMAX, NMAX )
      COMPLEX*16         WORKASAV( NMAX, NMAX )
      COMPLEX*16         WORKB( NMAX, MAXRHS )
      COMPLEX*16         WORKXACT( NMAX, MAXRHS )
      COMPLEX*16         WORKBSAV( NMAX, MAXRHS )
      COMPLEX*16         WORKX( NMAX, MAXRHS )
      COMPLEX*16         WORKAFAC( NMAX, NMAX )
      COMPLEX*16         WORKAINV( NMAX, NMAX )
      COMPLEX*16         WORKARF( (NMAX*(NMAX+1))/2 )
      COMPLEX*16         WORKAP( (NMAX*(NMAX+1))/2 )
      COMPLEX*16         WORKARFINV( (NMAX*(NMAX+1))/2 )
      COMPLEX*16         Z_WORK_ZLATMS( 3 * NMAX )
      COMPLEX*16         Z_WORK_ZPOT02( NMAX, MAXRHS )
      COMPLEX*16         Z_WORK_ZPOT03( NMAX, NMAX )
      double             D_WORK_ZLATMS( NMAX );
      double             D_WORK_ZLANHE( NMAX );
      double             D_WORK_ZPOT01( NMAX );
      double             D_WORK_ZPOT02( NMAX );
      double             D_WORK_ZPOT03( NMAX );
      // ..
      // .. External Functions ..
      double             DLAMCH, DSECND;
      // EXTERNAL DLAMCH, DSECND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ILAVER, ZDRVRFP, ZDRVRF1, ZDRVRF2, ZDRVRF3, ZDRVRF4
      // ..
      // .. Executable Statements ..

      S1 = DSECND( )
      FATAL = .FALSE.

      // Read a dummy line.

      READ( NIN, FMT = * )

      // Report LAPACK version tag (e.g. LAPACK-3.2.0)

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH

      // Read the values of N

      READ( NIN, FMT = * )NN
      if ( NN.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      } else if ( NN.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      for (I = 1; I <= NN; I++) { // 10
         if ( NVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' M  ', NVAL( I ), 0
            FATAL = .TRUE.
         } else if ( NVAL( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9995 )' M  ', NVAL( I ), NMAX
            FATAL = .TRUE.
         }
   10 CONTINUE
      IF( NN.GT.0 ) WRITE( NOUT, FMT = 9993 )'N   ', ( NVAL( I ), I = 1, NN )

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS
      if ( NNS.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1
         NNS = 0
         FATAL = .TRUE.
      } else if ( NNS.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN
         NNS = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS )
      for (I = 1; I <= NNS; I++) { // 30
         if ( NSVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0
            FATAL = .TRUE.
         } else if ( NSVAL( I ).GT.MAXRHS ) {
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS
            FATAL = .TRUE.
         }
   30 CONTINUE
      IF( NNS.GT.0 ) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )

      // Read the matrix types

      READ( NIN, FMT = * )NNT
      if ( NNT.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NMA', NNT, 1
         NNT = 0
         FATAL = .TRUE.
      } else if ( NNT.GT.NTYPES ) {
         WRITE( NOUT, FMT = 9995 )' NMA', NNT, NTYPES
         NNT = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( NTVAL( I ), I = 1, NNT )
      for (I = 1; I <= NNT; I++) { // 320
         if ( NTVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )'TYPE', NTVAL( I ), 0
            FATAL = .TRUE.
         } else if ( NTVAL( I ).GT.NTYPES ) {
            WRITE( NOUT, FMT = 9995 )'TYPE', NTVAL( I ), NTYPES
            FATAL = .TRUE.
         }
  320 CONTINUE
      IF( NNT.GT.0 ) WRITE( NOUT, FMT = 9993 )'TYPE', ( NTVAL( I ), I = 1, NNT )

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9992 )THRESH

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 )
         STOP
      }

      // Calculate and print the machine dependent constants.

      EPS = DLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9991 )'underflow', EPS
      EPS = DLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9991 )'overflow ', EPS
      EPS = DLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9991 )'precision', EPS
      WRITE( NOUT, FMT = * )

      // Test the error exit of:

      IF( TSTERR ) CALL ZERRRFP( NOUT )

*    Test the routines: zpftrf, zpftri, zpftrs (as in ZDRVPO).
*    This also tests the routines: ztfsm, ztftri, ztfttr, ztrttf.

      zdrvrfp(NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, WORKA, WORKASAV, WORKAFAC, WORKAINV, WORKB, WORKBSAV, WORKXACT, WORKX, WORKARF, WORKARFINV, Z_WORK_ZLATMS, Z_WORK_ZPOT02, Z_WORK_ZPOT03, D_WORK_ZLATMS, D_WORK_ZLANHE, D_WORK_ZPOT01, D_WORK_ZPOT02, D_WORK_ZPOT03 );

*    Test the routine: zlanhf

      zdrvrf1(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, D_WORK_ZLANHE );

*    Test the conversion routines:
        // zhfttp, ztpthf, ztfttr, ztrttf, ztrttp and ztpttr.

      zdrvrf2(NOUT, NN, NVAL, WORKA, NMAX, WORKARF, WORKAP, WORKASAV );

*    Test the routine: ztfsm

      zdrvrf3(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, WORKAINV, WORKAFAC, D_WORK_ZLANHE, Z_WORK_ZPOT03, Z_WORK_ZPOT02 );


*    Test the routine: zhfrk

      zdrvrf4(NOUT, NN, NVAL, THRESH, WORKA, WORKAFAC, NMAX, WORKARF, WORKAINV, NMAX,D_WORK_ZLANHE);

      CLOSE ( NIN )
      S2 = DSECND( )
      WRITE( NOUT, FMT = 9998 )
      WRITE( NOUT, FMT = 9997 )S2 - S1

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9996 FORMAT( ' !! Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
 9995 FORMAT( ' !! Invalid input value: ', A4, '=', I6, '; must be <=', I6 )
 9994 FORMAT( /  ' Tests of the COMPLEX*16 LAPACK RFP routines ', / ' LAPACK VERSION ', I1, '.', I1, '.', I1, / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )

      // End of ZCHKRFP

      }
