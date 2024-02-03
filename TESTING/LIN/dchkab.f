      PROGRAM DCHKAB

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 132 ;
      int                MAXIN;
      const              MAXIN = 12 ;
      int                MAXRHS;
      const              MAXRHS = 16 ;
      int                MATMAX;
      const              MATMAX = 30 ;
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                LDAMAX;
      const              LDAMAX = NMAX ;
      // ..
      // .. Local Scalars ..
      bool               FATAL, TSTDRV, TSTERR;
      String             C1;
      String             C2;
      String             PATH;
      String             INTSTR;
      String             ALINE;
      int                I, IC, K, LDA, NM, NMATS, NNS, NRHS, NTYPES, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      double             EPS, S1, S2, THRESH;
      REAL               SEPS
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( MATMAX );
      int                IWORK( NMAX ), MVAL( MAXIN ), NSVAL( MAXIN );
      double             A( LDAMAX*NMAX, 2 ), B( NMAX*MAXRHS, 2 ), RWORK( NMAX ), WORK( NMAX*MAXRHS*2 );
      REAL               SWORK(NMAX*(NMAX+MAXRHS))
      // ..
      // .. External Functions ..
      double             DLAMCH, DSECND;
      bool               LSAME, LSAMEN;
      REAL               SLAMCH
      // EXTERNAL LSAME, LSAMEN, DLAMCH, DSECND, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, DDRVAB, DDRVAC, DERRAB, DERRAC, ILAVER
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               INTSTR / '0123456789' /
      // ..
      // .. Executable Statements ..

      S1 = DSECND( )
      LDA = NMAX
      FATAL = .FALSE.

      // Read a dummy line.

      READ( NIN, FMT = * )

      // Report values of parameters.

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH

      // Read the values of M

      READ( NIN, FMT = * )NM
      if ( NM.LT.1 ) {
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
         NM = 0
         FATAL = .TRUE.
      } else if ( NM.GT.MAXIN ) {
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
         NM = 0
         FATAL = .TRUE.
      }
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      for (I = 1; I <= NM; I++) { // 10
         if ( MVAL( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         } else if ( MVAL( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         }
      } // 10
      IF( NM.GT.0 ) WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM )

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
      } // 30
      IF( NNS.GT.0 ) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9992 )THRESH

      // Read the flag that indicates whether to test the driver routine.

      READ( NIN, FMT = * )TSTDRV

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      if ( FATAL ) {
         WRITE( NOUT, FMT = 9999 )
         STOP
      }

      // Calculate and print the machine dependent constants.

      SEPS = SLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9991 )'(single precision) underflow', SEPS
      SEPS = SLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9991 )'(single precision) overflow ', SEPS
      SEPS = SLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9991 )'(single precision) precision', SEPS
      WRITE( NOUT, FMT = * )

      EPS = DLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9991 )'(double          ) underflow', EPS;
      EPS = DLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9991 )'(double          ) overflow ', EPS;
      EPS = DLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9991 )'(double          ) precision', EPS;
      WRITE( NOUT, FMT = * )

      } // 80

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE
      PATH = ALINE( 1: 3 )
      NMATS = MATMAX
      I = 3
      } // 90
      I = I + 1
      if ( I.GT.72 ) {
         NMATS = MATMAX
         GO TO 130
      }
      IF( ALINE( I: I ).EQ.' ' ) GO TO 90
      NMATS = 0
      } // 100
      C1 = ALINE( I: I )
      for (K = 1; K <= 10; K++) { // 110
         if ( C1.EQ.INTSTR( K: K ) ) {
            IC = K - 1
            GO TO 120
         }
      } // 110
      GO TO 130
      } // 120
      NMATS = NMATS*10 + IC
      I = I + 1
      IF( I.GT.72 ) GO TO 130
      GO TO 100
      } // 130
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
      NRHS = NSVAL( 1 )

      // Check first character for correct precision.

      IF( .NOT.LSAME( C1, 'double          ' ) ) THEN;
         WRITE( NOUT, FMT = 9990 )PATH


      } else if ( NMATS.LE.0 ) {

         // Check for a positive number of tests requested.

         WRITE( NOUT, FMT = 9989 )PATH
         GO TO 140

      } else if ( LSAMEN( 2, C2, 'GE' ) ) {

         // GE:  general matrices

         NTYPES = 11
         alareq('DGE', NMATS, DOTYPE, NTYPES, NIN, NOUT );

         // Test the error exits

         IF( TSTERR ) CALL DERRAB( NOUT )

         if ( TSTDRV ) {
            ddrvab(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, IWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )'DSGESV'
         }

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // PO:  positive definite matrices

         NTYPES = 9
         alareq('DPO', NMATS, DOTYPE, NTYPES, NIN, NOUT );


         IF( TSTERR ) CALL DERRAC( NOUT )


         if ( TSTDRV ) {
            ddrvac(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, NOUT );
         } else {
            WRITE( NOUT, FMT = 9989 )PATH
         }
      } else {

      }

      // Go back to get another input line.

      GO TO 80

      // Branch to this line when the last record is read.

      } // 140
      CLOSE ( NIN )
      S2 = DSECND( )
      WRITE( NOUT, FMT = 9998 )
      WRITE( NOUT, FMT = 9997 )S2 - S1

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=', I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', I6 )
 9994 FORMAT( ' Tests of the double           LAPACK DSGESV/DSPOSV',; ' routines ', / ' LAPACK VERSION ', I1, '.', I1, '.', I1, / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ', 'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9990 FORMAT( / 1X, A6, ' routines were not tested' )
 9989 FORMAT( / 1X, A6, ' driver routines were not tested' )

      // End of DCHKAB

      }
