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

      CALL ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
      WRITE( NOUT, FMT = 9994 ) VERS_MAJOR, VERS_MINOR, VERS_PATCH

      // Read the values of M

      READ( NIN, FMT = * )NM
      IF( NM.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NM ', NM, 1
         NM = 0
         FATAL = .TRUE.
      ELSE IF( NM.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NM ', NM, MAXIN
         NM = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      DO 10 I = 1, NM
         IF( MVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( MVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9995 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   10 CONTINUE
      IF( NM.GT.0 ) WRITE( NOUT, FMT = 9993 )'M   ', ( MVAL( I ), I = 1, NM )

      // Read the values of NRHS

      READ( NIN, FMT = * )NNS
      IF( NNS.LT.1 ) THEN
         WRITE( NOUT, FMT = 9996 )' NNS', NNS, 1
         NNS = 0
         FATAL = .TRUE.
      ELSE IF( NNS.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9995 )' NNS', NNS, MAXIN
         NNS = 0
         FATAL = .TRUE.
      END IF
      READ( NIN, FMT = * )( NSVAL( I ), I = 1, NNS )
      DO 30 I = 1, NNS
         IF( NSVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9996 )'NRHS', NSVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NSVAL( I ).GT.MAXRHS ) THEN
            WRITE( NOUT, FMT = 9995 )'NRHS', NSVAL( I ), MAXRHS
            FATAL = .TRUE.
         END IF
   30 CONTINUE
      IF( NNS.GT.0 ) WRITE( NOUT, FMT = 9993 )'NRHS', ( NSVAL( I ), I = 1, NNS )

      // Read the threshold value for the test ratios.

      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9992 )THRESH

      // Read the flag that indicates whether to test the driver routine.

      READ( NIN, FMT = * )TSTDRV

      // Read the flag that indicates whether to test the error exits.

      READ( NIN, FMT = * )TSTERR

      IF( FATAL ) THEN
         WRITE( NOUT, FMT = 9999 )
         STOP
      END IF

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

   80 CONTINUE

      // Read a test path and the number of matrix types to use.

      READ( NIN, FMT = '(A72)', END = 140 )ALINE
      PATH = ALINE( 1: 3 )
      NMATS = MATMAX
      I = 3
   90 CONTINUE
      I = I + 1
      IF( I.GT.72 ) THEN
         NMATS = MATMAX
         GO TO 130
      END IF
      IF( ALINE( I: I ).EQ.' ' ) GO TO 90
      NMATS = 0
  100 CONTINUE
      C1 = ALINE( I: I )
      DO 110 K = 1, 10
         IF( C1.EQ.INTSTR( K: K ) ) THEN
            IC = K - 1
            GO TO 120
         END IF
  110 CONTINUE
      GO TO 130
  120 CONTINUE
      NMATS = NMATS*10 + IC
      I = I + 1
      IF( I.GT.72 ) GO TO 130
      GO TO 100
  130 CONTINUE
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
      NRHS = NSVAL( 1 )

      // Check first character for correct precision.

      IF( .NOT.LSAME( C1, 'double          ' ) ) THEN;
         WRITE( NOUT, FMT = 9990 )PATH


      ELSE IF( NMATS.LE.0 ) THEN

         // Check for a positive number of tests requested.

         WRITE( NOUT, FMT = 9989 )PATH
         GO TO 140

      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN

         // GE:  general matrices

         NTYPES = 11
         CALL ALAREQ( 'DGE', NMATS, DOTYPE, NTYPES, NIN, NOUT )

         // Test the error exits

         IF( TSTERR ) CALL DERRAB( NOUT )

         IF( TSTDRV ) THEN
            CALL DDRVAB( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, IWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )'DSGESV'
         END IF

      ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN

         // PO:  positive definite matrices

         NTYPES = 9
         CALL ALAREQ( 'DPO', NMATS, DOTYPE, NTYPES, NIN, NOUT )


         IF( TSTERR ) CALL DERRAC( NOUT )


         IF( TSTDRV ) THEN
            CALL DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, LDA, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ), WORK, RWORK, SWORK, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )PATH
         END IF
      ELSE

      END IF

      // Go back to get another input line.

      GO TO 80

      // Branch to this line when the last record is read.

  140 CONTINUE
      CLOSE ( NIN )
      S2 = DSECND( )
      WRITE( NOUT, FMT = 9998 )
      WRITE( NOUT, FMT = 9997 )S2 - S1

 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9998 FORMAT( / ' End of tests' )
 9997 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9996 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be >=',
     $      I6 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=',
     $      I6 )
 9994 FORMAT( ' Tests of the double           LAPACK DSGESV/DSPOSV',;
     $  ' routines ',
     $      / ' LAPACK VERSION ', I1, '.', I1, '.', I1,
     $      / / ' The following parameter values will be used:' )
 9993 FORMAT( 4X, A4, ':  ', 10I6, / 11X, 10I6 )
 9992 FORMAT( / ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / )
 9991 FORMAT( ' Relative machine ', A, ' is taken to be', D16.6 )
 9990 FORMAT( / 1X, A6, ' routines were not tested' )
 9989 FORMAT( / 1X, A6, ' driver routines were not tested' )

      // End of DCHKAB

      }
