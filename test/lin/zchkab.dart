      void main() {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

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
      bool               FATAL, TSTDRV, TSTERR;
      String             C1;
      String             C2;
      String             PATH;
      String             INTSTR;
      String             ALINE;
      int                I, IC, K, LDA, NM, NMATS, NNS, NRHS, NTYPES, VERS_MAJOR, VERS_MINOR, VERS_PATCH;
      double             EPS, S1, S2, THRESH;
      double               SEPS;
      bool               DOTYPE( MATMAX );
      int                IWORK( NMAX ), MVAL( MAXIN ), NSVAL( MAXIN );
      double             RWORK(NMAX);
      Complex         A( LDAMAX*NMAX, 2 ), B( NMAX*MAXRHS, 2 ), WORK( NMAX*MAXRHS*2 );
      Complex            SWORK(NMAX*(NMAX+MAXRHS));
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DSECND;
      //- bool               lsame, LSAMEN;
      //- REAL               SLAMCH;
      // EXTERNAL DLAMCH, DSECND, lsame, LSAMEN, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAREQ, ZDRVAB, ZDRVAC, ZERRAB, ZERRAC, ILAVER
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      // .. Data statements ..
      const INTSTR = '0123456789';

      S1 = DSECND( );
      LDA = NMAX;
      FATAL = false;

      // Read a dummy line.

      READ( NIN, FMT = * );

      // Report values of parameters.

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH );
      NOUT.println( 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH;

      // Read the values of M

      READ( NIN, FMT = * )NM;
      if ( NM < 1 ) {
         NOUT.println( 9996 )' NM ', NM, 1;
         NM = 0;
         FATAL = true;
      } else if ( NM > MAXIN ) {
         NOUT.println( 9995 )' NM ', NM, MAXIN;
         NM = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM );
      for (I = 1; I <= NM; I++) { // 10
         if ( MVAL( I ) < 0 ) {
            NOUT.println( 9996 )' M  ', MVAL( I ), 0;
            FATAL = true;
         } else if ( MVAL( I ) > NMAX ) {
            NOUT.println( 9995 )' M  ', MVAL( I ), NMAX;
            FATAL = true;
         }
      } // 10
      if (NM > 0) NOUT.println( 9993 )'M   ', ( MVAL( I ), I = 1, NM );

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS;
      if ( NNS < 1 ) {
         NOUT.println( 9996 )' NNS', NNS, 1;
         NNS = 0;
         FATAL = true;
      } else if ( NNS > MAXIN ) {
         NOUT.println( 9995 )' NNS', NNS, MAXIN;
         NNS = 0;
         FATAL = true;
      }
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS );
      for (I = 1; I <= NNS; I++) { // 30
         if ( NSVAL( I ) < 0 ) {
            NOUT.println( 9996 )'NRHS', NSVAL( I ), 0;
            FATAL = true;
         } else if ( NSVAL( I ) > MAXRHS ) {
            NOUT.println( 9995 )'NRHS', NSVAL( I ), MAXRHS;
            FATAL = true;
         }
      } // 30
      if (NNS > 0) NOUT.println( 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS );

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH;
      NOUT.println( 9992 )THRESH;

      // Read the flag that indicates whether to test the driver routine.

      READ( NIN, FMT = * )TSTDRV;

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR;

      if ( FATAL ) {
         NOUT.println( 9999 );
         STOP;
      }

      // Calculate and print the machine dependent constants.

      SEPS = SLAMCH( 'Underflow threshold' );
      NOUT.println( 9991 )'(single precision) underflow', SEPS;
      SEPS = SLAMCH( 'Overflow threshold' );
      NOUT.println( 9991 )'(single precision) overflow ', SEPS;
      SEPS = SLAMCH( 'Epsilon' );
      NOUT.println( 9991 )'(single precision) precision', SEPS;
      NOUT.println( * );

      EPS = dlamch( 'Underflow threshold' );
      NOUT.println( 9991 )'(double          ) underflow', EPS;
      EPS = dlamch( 'Overflow threshold' );
      NOUT.println( 9991 )'(double          ) overflow ', EPS;
      EPS = dlamch( 'Epsilon' );
      NOUT.println( 9991 )'(double          ) precision', EPS;
      NOUT.println( * );

      } // 80

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE;
      PATH = ALINE( 1: 3 );
      NMATS = MATMAX;
      I = 3;
      } // 90
      I++;
      if ( I > 72 ) {
         NMATS = MATMAX;
         GO TO 130;
      }
      if( ALINE( I: I ) == ' ' ) GO TO 90;
      NMATS = 0;
      } // 100
      C1 = ALINE( I: I );
      for (K = 1; K <= 10; K++) { // 110
         if ( C1 == INTSTR( K: K ) ) {
            IC = K - 1;
            GO TO 120;
         }
      } // 110
      GO TO 130;
      } // 120
      NMATS = NMATS*10 + IC;
      I++;
      if (I > 72) GO TO 130;
      GO TO 100;
      } // 130
      C1 = PATH( 1: 1 );
      C2 = PATH.substring( 1, 3 );
      final NRHS = NSVAL[1];
      final NRHS = NSVAL[1];

      // Check first character for correct precision.

      if ( !lsame( C1, 'Zomplex precision' ) ) {
            NOUT.println( 9990 )PATH;

      } else if ( NMATS <= 0 ) {

         // Check for a positive number of tests requested.

         NOUT.println( 9990 )'ZCGESV';
         GO TO 140;

      } else if ( lsamen( 2, C2, 'GE' ) ) {

         // GE:  general matrices

      NTYPES = 11;
      alareq('ZGE', NMATS, DOTYPE, NTYPES, NIN, NOUT );

         // Test the error exits

         if (TSTERR) zerrab( NOUT );

         if ( TSTDRV ) {
            zdrvab(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, IWORK, NOUT );
         } else {
            NOUT.println( 9989 )'ZCGESV';
         }

      } else if ( lsamen( 2, C2, 'PO' ) ) {

         // PO:  positive definite matrices

         NTYPES = 9;
         alareq('DPO', NMATS, DOTYPE, NTYPES, NIN, NOUT );

         if (TSTERR) zerrac( NOUT );


         if ( TSTDRV ) {
            zdrvac(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, NOUT );
         } else {
            NOUT.println( 9989 )'ZCPOSV';
         }

      } else {

      }

      // Go back to get another input line.

      GO TO 80;

      // Branch to this line when the last record is read.

      } // 140
      CLOSE ( NIN );
      S2 = DSECND( );
      NOUT.println( 9998 );
      NOUT.println( 9997 )S2 - S1;

 9999 FORMAT('\n Execution not attempted due to input errors' );
 9998 FORMAT('\n End of tests' );
 9997 FORMAT( ' Total time used = ${.f12_2} seconds\n');
 9996 FORMAT( ' Invalid input value: ${.a4}=${.i6}; must be >=', I6 )
 9995 FORMAT( ' Invalid input value: ${.a4}=${.i6}; must be <=', I6 )
 9994 FORMAT( ' Tests of the Complex LAPACK ZCGESV/ZCPOSV routines \n LAPACK VERSION ${.i1}.${.i1}.', I1, / / ' The following parameter values will be used:' );
 9993 FORMAT('    ${.a4}:  ', 10I6, / 11X, 10I6 );
 9992 FORMAT('\n Routines pass computational tests if test ratio is less than', F8.2, / );
 9991 FORMAT( ' Relative machine ${} is taken to be', D16.6 );
 9990 FORMAT('\n ${.a6} routines were not tested' );
 9989 FORMAT('\n ${.a6} driver routines were not tested' );
      }
