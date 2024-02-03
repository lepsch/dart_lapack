      void main() {
// Test program for the Complex          Level 2 Blas.

// The program must be driven by a short data file. The first 17 records
// of the file are read using list-directed input, the last 17 records
// are read using the format ( A12, L2 ). An annotated example of a data
// file can be obtained by deleting the first 3 characters from the
// following 34 lines:
// 'CBLAT2.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
// -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF < 0)
// F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
// F        LOGICAL FLAG, T TO STOP ON FAILURES.
// T        LOGICAL FLAG, T TO TEST ERROR EXITS.
// 2        0 TO TEST COLUMN-MAJOR, 1 TO TEST ROW-MAJOR, 2 TO TEST BOTH
// 16.0     THRESHOLD VALUE OF TEST RATIO
// 6                 NUMBER OF VALUES OF N
// 0 1 2 3 5 9       VALUES OF N
// 4                 NUMBER OF VALUES OF K
// 0 1 2 4           VALUES OF K
// 4                 NUMBER OF VALUES OF INCX AND INCY
// 1 2 -1 -2         VALUES OF INCX AND INCY
// 3                 NUMBER OF VALUES OF ALPHA
// (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
// 3                 NUMBER OF VALUES OF BETA
// (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
// cblas_zgemv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zgbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zhemv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zhbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zhpmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztrmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztpmv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztrsv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztbsv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ztpsv  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zgerc  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zgeru  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zher   T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zhpr   T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zher2  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_zhpr2  T PUT F FOR NO TEST. SAME COLUMNS.

      // See:

         // Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
         // An  extended  set of Fortran  Basic Linear Algebra Subprograms.

         // Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
         // and  Computer Science  Division,  Argonne  National Laboratory,
         // 9700 South Cass Avenue, Argonne, Illinois 60439, US.

         // Or

         // NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
         // Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
         // OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
         // Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.


// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                NSUBS;
      const              NSUBS = 17 ;
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO, RHALF, RONE;
      const              RZERO = 0.0, RHALF = 0.5, RONE = 1.0 ;
      int                NMAX, INCMAX;
      const              NMAX = 65, INCMAX = 2 ;
      int                NINMAX, NIDMAX, NKBMAX, NALMAX, NBEMAX;
      const              NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      double             EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB, NTRA, LAYOUT;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR, CORDER, RORDER;
      String             TRANS;
      String             SNAMET;
      String             SNAPS;
      // .. Local Arrays ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BET( NBEMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( 2*NMAX );
      double             G( NMAX );
      int                IDIM( NIDMAX ), INC( NINMAX ), KB( NKBMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      //- double             DDIFF;
      //- bool               LZE;
      // EXTERNAL DDIFF, LZE
      // .. External Subroutines ..
      // EXTERNAL ZCHK1, ZCHK2, ZCHK3, ZCHK4, ZCHK5, ZCHK6, CZ2CHKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               OK;
      String             SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // COMMON /SRNAMC/SRNAMT
      // .. Data statements ..
      const SNAMES = ['cblas_zgemv ', 'cblas_zgbmv ', 'cblas_zhemv ','cblas_zhbmv ','cblas_zhpmv ', 'cblas_ztrmv ','cblas_ztbmv ','cblas_ztpmv ', 'cblas_ztrsv ','cblas_ztbsv ','cblas_ztpsv ', 'cblas_zgerc ','cblas_zgeru ','cblas_zher  ', 'cblas_zhpr  ','cblas_zher2 ','cblas_zhpr2 '];
      // .. Executable Statements ..

      NOUTC = NOUT;

      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SNAPS;
      READ( NIN, FMT = * )NTRA;
      TRACE = NTRA >= 0;
      if ( TRACE ) {
         OPEN( NTRA, FILE = SNAPS );
      }
      // Read the flag that directs rewinding of the snapshot file.
      READ( NIN, FMT = * )REWI;
      REWI = REWI && TRACE;
      // Read the flag that directs stopping on any failure.
      READ( NIN, FMT = * )SFATAL;
      // Read the flag that indicates whether error exits are to be tested.
      READ( NIN, FMT = * )TSTERR;
      // Read the flag that indicates whether row-major data layout to be tested.
      READ( NIN, FMT = * )LAYOUT;
      // Read the threshold value of the test ratio
      READ( NIN, FMT = * )THRESH;

      // Read and check the parameter values for the tests.

      // Values of N
      READ( NIN, FMT = * )NIDIM;
      if ( NIDIM < 1 || NIDIM > NIDMAX ) {
         WRITE( NOUT, FMT = 9997 )'N', NIDMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM );
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ) < 0 || IDIM( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX;
            GO TO 230;
         }
      } // 10
      // Values of K
      READ( NIN, FMT = * )NKB;
      if ( NKB < 1 || NKB > NKBMAX ) {
         WRITE( NOUT, FMT = 9997 )'K', NKBMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( KB( I ), I = 1, NKB );
      for (I = 1; I <= NKB; I++) { // 20
         if ( KB( I ) < 0 ) {
            WRITE( NOUT, FMT = 9995 );
            GO TO 230;
         }
      } // 20
      // Values of INCX and INCY
      READ( NIN, FMT = * )NINC;
      if ( NINC < 1 || NINC > NINMAX ) {
         WRITE( NOUT, FMT = 9997 )'INCX AND INCY', NINMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( INC( I ), I = 1, NINC );
      for (I = 1; I <= NINC; I++) { // 30
         if ( INC( I ) == 0 || ( INC( I ) ).abs() > INCMAX ) {
            WRITE( NOUT, FMT = 9994 )INCMAX;
            GO TO 230;
         }
      } // 30
      // Values of ALPHA
      READ( NIN, FMT = * )NALF;
      if ( NALF < 1 || NALF > NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF );
      // Values of BETA
      READ( NIN, FMT = * )NBET;
      if ( NBET < 1 || NBET > NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX;
         GO TO 230;
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET );

      // Report values of parameters.

      WRITE( NOUT, FMT = 9993 );
      WRITE( NOUT, FMT = 9992 )( IDIM( I ), I = 1, NIDIM );
      WRITE( NOUT, FMT = 9991 )( KB( I ), I = 1, NKB );
      WRITE( NOUT, FMT = 9990 )( INC( I ), I = 1, NINC );
      WRITE( NOUT, FMT = 9989 )( ALF( I ), I = 1, NALF );
      WRITE( NOUT, FMT = 9988 )( BET( I ), I = 1, NBET );
      if ( !TSTERR ) {
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = 9980 );
      }
      WRITE( NOUT, FMT = * );
      WRITE( NOUT, FMT = 9999 )THRESH;
      WRITE( NOUT, FMT = * );
      RORDER = false;
      CORDER = false;
      if (LAYOUT == 2) {
         RORDER = true;
         CORDER = true;
         WRITE( *, FMT = 10002 );
      } else if (LAYOUT == 1) {
         RORDER = true;
         WRITE( *, FMT = 10001 );
      } else if (LAYOUT == 0) {
         CORDER = true;
         WRITE( *, FMT = 10000 );
      }
      WRITE( *, FMT = * );

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      for (I = 1; I <= NSUBS; I++) { // 40
         LTEST( I ) = false;
      } // 40
   50 READ( NIN, FMT = 9984, END = 80 )SNAMET, LTESTT;
      for (I = 1; I <= NSUBS; I++) { // 60
         if( SNAMET == SNAMES( I ) ) GO TO 70;
      } // 60
      WRITE( NOUT, FMT = 9986 )SNAMET;
      STOP;
   70 LTEST( I ) = LTESTT;
      GO TO 50;

      } // 80
      CLOSE ( NIN );

      // Compute EPS (the machine precision).

      EPS = RONE;
      } // 90
      if( DDIFF( RONE + EPS, RONE ) == RZERO ) GO TO 100;
      EPS = RHALF*EPS;
      GO TO 90;
      } // 100
      EPS = EPS + EPS;
      WRITE( NOUT, FMT = 9998 )EPS;

      // Check the reliability of ZMVCH using exact data.

      N = min( 32, NMAX );
      for (J = 1; J <= N; J++) { // 120
         for (I = 1; I <= N; I++) { // 110
            A( I, J ) = max( I - J + 1, 0 );
         } // 110
         X( J ) = J;
         Y( J ) = ZERO;
      } // 120
      for (J = 1; J <= N; J++) { // 130
         YY( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
      } // 130
      // YY holds the exact result. On exit from CMVCH YT holds
      // the result computed by CMVCH.
      TRANS = 'N';
      zmvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LZE( YY, YT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR;
         STOP;
      }
      TRANS = 'T';
      zmvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LZE( YY, YT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR;
         STOP;
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 210
         WRITE( NOUT, FMT = * );
         if ( !LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM );
         } else {
            SRNAMT = SNAMES( ISNUM );
            // Test error exits.
            if ( TSTERR ) {
               cz2chke(SNAMES( ISNUM ) );
               WRITE( NOUT, FMT = * );
            }
            // Test computations.
            INFOT = 0;
            OK = true;
            FATAL = false;
            GO TO ( 140, 140, 150, 150, 150, 160, 160, 160, 160, 160, 160, 170, 170, 180, 180, 190, 190 )ISNUM;
            // Test ZGEMV, 01, and ZGBMV, 02.
  140       IF (CORDER) THEN;
            zchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 0 );
            }
            if (RORDER) {
            zchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 1 );
            }
            GO TO 200;
            // Test ZHEMV, 03, ZHBMV, 04, and ZHPMV, 05.
  150      IF (CORDER) THEN;
           zchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 0 );
           }
           if (RORDER) {
           zchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 1 );
           }
            GO TO 200;
            // Test ZTRMV, 06, ZTBMV, 07, ZTPMV, 08,
            // ZTRSV, 09, ZTBSV, 10, and ZTPSV, 11.
  160      IF (CORDER) THEN;
           zchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200;
            // Test ZGERC, 12, ZGERU, 13.
  170      IF (CORDER) THEN;
           zchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200;
            // Test ZHER, 14, and ZHPR, 15.
  180      IF (CORDER) THEN;
           zchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200;
            // Test ZHER2, 16, and ZHPR2, 17.
  190      IF (CORDER) THEN;
           zchk6(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk6(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }

  200       IF( FATAL && SFATAL ) GO TO 220;
         }
      } // 210
      WRITE( NOUT, FMT = 9982 );
      GO TO 240;

      } // 220
      WRITE( NOUT, FMT = 9981 );
      GO TO 240;

      } // 230
      WRITE( NOUT, FMT = 9987 );

      } // 240
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT );
      STOP;

10002 FORMAT( ' COLUMN-MAJOR AND ROW-MAJOR DATA LAYOUTS ARE TESTED' )
10001 FORMAT( ' ROW-MAJOR DATA LAYOUT IS TESTED' )
10000 FORMAT( ' COLUMN-MAJOR DATA LAYOUT IS TESTED' )
 9999 FORMAT(' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 );
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 );
 9997 FORMAT(' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 );
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 );
 9995 FORMAT( ' VALUE OF K IS LESS THAN 0' );
 9994 FORMAT( ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ', I2 );
 9993 FORMAT(' TESTS OF THE Complex      LEVEL 2 BLAS', //' THE F', 'OLLOWING PARAMETER VALUES WILL BE USED:' );
 9992 FORMAT( '   FOR N              ', 9I6 );
 9991 FORMAT( '   FOR K              ', 7I6 );
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 );
 9989 FORMAT( '   FOR ALPHA          ', 7('(', F4.1, ',', F4.1, ')  ', : ) );
 9988 FORMAT( '   FOR BETA           ', 7('(', F4.1, ',', F4.1, ')  ', : ) );
 9987 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' );
 9986 FORMAT(' SUBPROGRAM NAME ',A12, ' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' );
 9985 FORMAT(' ERROR IN CMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' CMVCH WAS CALLED WITH TRANS = ', A1, ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '.', / ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.' , /' ******* TESTS ABANDONED *******' );
 9984 FORMAT( A12, L2 );
 9983 FORMAT( 1X,A12, ' WAS NOT TESTED' );
 9982 FORMAT( /' END OF TESTS' );
 9981 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' );
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' );

      // End of ZBLAT2.

      }
      void zchk1(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, IORDER ) {

// Tests CGEMV and CGBMV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex        ZERO, HALF;
      const             ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      double            RZERO;
      const             RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, BETA, BLS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA, LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK, NL, NS;
      bool               BANDED, FULL, NULL, RESET, SAME, TRAN;
      String             TRANS, TRANSS;
      String             CTRANS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZGBMV, CZGEMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      const ICH = 'NTC';
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ) == 'e';
      BANDED = SNAME( 9: 9 ) == 'b';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 11;
      } else if ( BANDED ) {
         NARGS = 13;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN );
         ND = N/2 + 1;

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = max( N - ND, 0 );
            IF( IM == 2 ) M = min( N + ND, NMAX );

            if ( BANDED ) {
               NK = NKB;
            } else {
               NK = 1;
            }
            for (IKU = 1; IKU <= NK; IKU++) { // 100
               if ( BANDED ) {
                  KU = KB( IKU );
                  KL = max( KU - 1, 0 );
               } else {
                  KU = N - 1;
                  KL = M - 1;
               }
               // Set LDA to 1 more than minimum value if room.
               if ( BANDED ) {
                  LDA = KL + KU + 1;
               } else {
                  LDA = M;
               }
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 100;
               LAA = LDA*N;
               NULL = N <= 0 || M <= 0;

               // Generate the matrix A.

               TRANSL = ZERO;
               zmake(SNAME( 8: 9 ), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL );

               for (IC = 1; IC <= 3; IC++) { // 90
                  TRANS = ICH( IC: IC );
                  if (TRANS == 'N') {
                     CTRANS = '  CblasNoTrans';
                  } else if (TRANS == 'T') {
                     CTRANS = '    CblasTrans';
                  } else {
                     CTRANS = 'CblasConjTrans';
                  }
                  TRAN = TRANS == 'T' || TRANS == 'C';

                  if ( TRAN ) {
                     ML = N;
                     NL = M;
                  } else {
                     ML = M;
                     NL = N;
                  }

                  for (IX = 1; IX <= NINC; IX++) { // 80
                     INCX = INC( IX );
                     LX = ( INCX ).abs()*NL;

                     // Generate the vector X.

                     TRANSL = HALF;
                     zmake('ge', ' ', ' ', 1, NL, X, 1, XX, ( INCX ).abs(), 0, NL - 1, RESET, TRANSL );
                     if ( NL > 1 ) {
                        X( NL/2 ) = ZERO;
                        XX( 1 + ( INCX ).abs()*( NL/2 - 1 ) ) = ZERO;
                     }

                     for (IY = 1; IY <= NINC; IY++) { // 70
                        INCY = INC( IY );
                        LY = ( INCY ).abs()*ML;

                        for (IA = 1; IA <= NALF; IA++) { // 60
                           ALPHA = ALF( IA );

                           for (IB = 1; IB <= NBET; IB++) { // 50
                              BETA = BET( IB );

                              // Generate the vector Y.

                              TRANSL = ZERO;
                              zmake('ge', ' ', ' ', 1, ML, Y, 1, YY, ( INCY ).abs(), 0, ML - 1, RESET, TRANSL );

                              NC = NC + 1;

                              // Save every datum before calling the
                              // subroutine.

                              TRANSS = TRANS;
                              MS = M;
                              NS = N;
                              KLS = KL;
                              KUS = KU;
                              ALS = ALPHA;
                              for (I = 1; I <= LAA; I++) { // 10
                                 AS( I ) = AA( I );
                              } // 10
                              LDAS = LDA;
                              for (I = 1; I <= LX; I++) { // 20
                                 XS( I ) = XX( I );
                              } // 20
                              INCXS = INCX;
                              BLS = BETA;
                              for (I = 1; I <= LY; I++) { // 30
                                 YS( I ) = YY( I );
                              } // 30
                              INCYS = INCY;

                              // Call the subroutine.

                              if ( FULL ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CTRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 czgemv(IORDER, TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              } else if ( BANDED ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, CTRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 czgbmv(IORDER, TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              }

                             // Check if error-exit was taken incorrectly.

                              if ( !OK ) {
                                 WRITE( NOUT, FMT = 9993 );
                                 FATAL = true;
                                 GO TO 130;
                              }

                              // See what data changed inside subroutines.

         // IF(TRANS != 'C' || (INCX > 0 && INCY > 0)) THEN
                              ISAME( 1 ) = TRANS == TRANSS;
                              ISAME( 2 ) = MS == M;
                              ISAME( 3 ) = NS == N;
                              if ( FULL ) {
                                 ISAME( 4 ) = ALS == ALPHA;
                                 ISAME( 5 ) = LZE( AS, AA, LAA );
                                 ISAME( 6 ) = LDAS == LDA;
                                 ISAME( 7 ) = LZE( XS, XX, LX );
                                 ISAME( 8 ) = INCXS == INCX;
                                 ISAME( 9 ) = BLS == BETA;
                                 if ( NULL ) {
                                    ISAME( 10 ) = LZE( YS, YY, LY );
                                 } else {
                                    ISAME( 10 ) = LZERES( 'ge', ' ', 1, ML, YS, YY, ( INCY ).abs() );
                                 }
                                 ISAME( 11 ) = INCYS == INCY;
                              } else if ( BANDED ) {
                                 ISAME( 4 ) = KLS == KL;
                                 ISAME( 5 ) = KUS == KU;
                                 ISAME( 6 ) = ALS == ALPHA;
                                 ISAME( 7 ) = LZE( AS, AA, LAA );
                                 ISAME( 8 ) = LDAS == LDA;
                                 ISAME( 9 ) = LZE( XS, XX, LX );
                                 ISAME( 10 ) = INCXS == INCX;
                                 ISAME( 11 ) = BLS == BETA;
                                 if ( NULL ) {
                                    ISAME( 12 ) = LZE( YS, YY, LY );
                                 } else {
                                    ISAME( 12 ) = LZERES( 'ge', ' ', 1, ML, YS, YY, ( INCY ).abs() );
                                 }
                                 ISAME( 13 ) = INCYS == INCY;
                              }

                              // If data was incorrectly changed, report
                              // and return.

                              SAME = true;
                              for (I = 1; I <= NARGS; I++) { // 40
                                 SAME = SAME && ISAME( I );
                                 if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                              } // 40
                              if ( !SAME ) {
                                 FATAL = true;
                                 GO TO 130;
                              }

                              if ( !NULL ) {

                                 // Check the result.

                                 zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                                 ERRMAX = max( ERRMAX, ERR );
                                 // If got really bad answer, report and
                                 // return.
                                 if (FATAL) GO TO 130;
                              } else {
                                 // Avoid repeating tests with M <= 0 or
                                 // N <= 0.
                                 GO TO 110;
                              }
                           // END IF

                           } // 50

                        } // 60

                     } // 70

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 140;

      } // 130
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CTRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CTRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
      }

      } // 140
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 4( I3, ',' ), '(', F4.1, ',', F4.1, '), A,',/ 10x, I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' );
 9994 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,',/ 10x, I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' );
 9993 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of ZCHK1.

      }
      void zchk2(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, IORDER ) {

// Tests CHEMV, CHBMV and CHPMV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, BETA, BLS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IK, IN, INCX, INCXS, INCY, INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHBMV, CZHEMV, CZHPMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ) == 'e';
      BANDED = SNAME( 9: 9 ) == 'b';
      PACKED = SNAME( 9: 9 ) == 'p';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 10;
      } else if ( BANDED ) {
         NARGS = 11;
      } else if ( PACKED ) {
         NARGS = 9;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN );

         if ( BANDED ) {
            NK = NKB;
         } else {
            NK = 1;
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK );
            } else {
               K = N - 1;
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1;
            } else {
               LDA = N;
            }
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2;
            } else {
               LAA = LDA*N;
            }
            NULL = N <= 0;

            for (IC = 1; IC <= 2; IC++) { // 90
               UPLO = ICH( IC: IC );
               if (UPLO == 'U') {
                  CUPLO = '    CblasUpper';
               } else {
                  CUPLO = '    CblasLower';
               }

               // Generate the matrix A.

               TRANSL = ZERO;
               zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

               for (IX = 1; IX <= NINC; IX++) { // 80
                  INCX = INC( IX );
                  LX = ( INCX ).abs()*N;

                  // Generate the vector X.

                  TRANSL = HALF;
                  zmake('ge', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     X( N/2 ) = ZERO;
                     XX( 1 + ( INCX ).abs()*( N/2 - 1 ) ) = ZERO;
                  }

                  for (IY = 1; IY <= NINC; IY++) { // 70
                     INCY = INC( IY );
                     LY = ( INCY ).abs()*N;

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA );

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB );

                           // Generate the vector Y.

                           TRANSL = ZERO;
                           zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );

                           NC = NC + 1;

                           // Save every datum before calling the
                           // subroutine.

                           UPLOS = UPLO;
                           NS = N;
                           KS = K;
                           ALS = ALPHA;
                           for (I = 1; I <= LAA; I++) { // 10
                              AS( I ) = AA( I );
                           } // 10
                           LDAS = LDA;
                           for (I = 1; I <= LX; I++) { // 20
                              XS( I ) = XX( I );
                           } // 20
                           INCXS = INCX;
                           BLS = BETA;
                           for (I = 1; I <= LY; I++) { // 30
                              YS( I ) = YY( I );
                           } // 30
                           INCYS = INCY;

                           // Call the subroutine.

                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, LDA, INCX, BETA, INCY;
                              IF( REWI ) REWIND NTRA;
                              czhemv(IORDER, UPLO, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, K, ALPHA, LDA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              czhbmv(IORDER, UPLO, N, K, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, N, ALPHA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              czhpmv(IORDER, UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( !OK ) {
                              WRITE( NOUT, FMT = 9992 );
                              FATAL = true;
                              GO TO 120;
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = UPLO == UPLOS;
                           ISAME( 2 ) = NS == N;
                           if ( FULL ) {
                              ISAME( 3 ) = ALS == ALPHA;
                              ISAME( 4 ) = LZE( AS, AA, LAA );
                              ISAME( 5 ) = LDAS == LDA;
                              ISAME( 6 ) = LZE( XS, XX, LX );
                              ISAME( 7 ) = INCXS == INCX;
                              ISAME( 8 ) = BLS == BETA;
                              if ( NULL ) {
                                 ISAME( 9 ) = LZE( YS, YY, LY );
                              } else {
                                 ISAME( 9 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME( 10 ) = INCYS == INCY;
                           } else if ( BANDED ) {
                              ISAME( 3 ) = KS == K;
                              ISAME( 4 ) = ALS == ALPHA;
                              ISAME( 5 ) = LZE( AS, AA, LAA );
                              ISAME( 6 ) = LDAS == LDA;
                              ISAME( 7 ) = LZE( XS, XX, LX );
                              ISAME( 8 ) = INCXS == INCX;
                              ISAME( 9 ) = BLS == BETA;
                              if ( NULL ) {
                                 ISAME( 10 ) = LZE( YS, YY, LY );
                              } else {
                                 ISAME( 10 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME( 11 ) = INCYS == INCY;
                           } else if ( PACKED ) {
                              ISAME( 3 ) = ALS == ALPHA;
                              ISAME( 4 ) = LZE( AS, AA, LAA );
                              ISAME( 5 ) = LZE( XS, XX, LX );
                              ISAME( 6 ) = INCXS == INCX;
                              ISAME( 7 ) = BLS == BETA;
                              if ( NULL ) {
                                 ISAME( 8 ) = LZE( YS, YY, LY );
                              } else {
                                 ISAME( 8 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ( INCY ).abs() );
                              }
                              ISAME( 9 ) = INCYS == INCY;
                           }

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 40
                              SAME = SAME && ISAME( I );
                              if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                           } // 40
                           if ( !SAME ) {
                              FATAL = true;
                              GO TO 120;
                           }

                           if ( !NULL ) {

                              // Check the result.

                              zmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                              ERRMAX = max( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 120;
                           } else {
                              // Avoid repeating tests with N <= 0
                              GO TO 110;
                           }

                        } // 50

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, K, ALPHA, LDA, INCX, BETA, INCY;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CUPLO, N, ALPHA, INCX, BETA, INCY;
      }

      } // 130
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), AP, X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' );
 9994 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' );
 9993 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), A,', I3, ', X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), ', 'Y,', I2, ') .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CZHK2.

      }
      void zchk3(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, XT, G, Z, IORDER ) {

// Tests ZTRMV, ZTBMV, ZTPMV, ZTRSV, ZTBSV and ZTPSV.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XT( NMAX ), XX( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      Complex         TRANSL;
      double             ERR, ERRMAX;
      int                I, ICD, ICT, ICU, IK, IN, INCX, INCXS, IX, K, KS, LAA, LDA, LDAS, LX, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             DIAG, DIAGS, TRANS, TRANSS, UPLO, UPLOS;
      String             CUPLO,CTRANS,CDIAG;
      String             ICHD, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZMAKE, ZMVCH, CZTBMV, CZTBSV, CZTPMV, CZTPSV, CZTRMV, CZTRSV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN';
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ) == 'r';
      BANDED = SNAME( 9: 9 ) == 'b';
      PACKED = SNAME( 9: 9 ) == 'p';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 8;
      } else if ( BANDED ) {
         NARGS = 9;
      } else if ( PACKED ) {
         NARGS = 7;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;
      // Set up zero vector for ZMVCH.
      for (I = 1; I <= NMAX; I++) { // 10
         Z( I ) = ZERO;
      } // 10

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN );

         if ( BANDED ) {
            NK = NKB;
         } else {
            NK = 1;
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK );
            } else {
               K = N - 1;
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1;
            } else {
               LDA = N;
            }
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2;
            } else {
               LAA = LDA*N;
            }
            NULL = N <= 0;

            for (ICU = 1; ICU <= 2; ICU++) { // 90
               UPLO = ICHU( ICU: ICU );
               if (UPLO == 'U') {
                  CUPLO = '    CblasUpper';
               } else {
                  CUPLO = '    CblasLower';
               }

               for (ICT = 1; ICT <= 3; ICT++) { // 80
                  TRANS = ICHT( ICT: ICT );
                  if (TRANS == 'N') {
                     CTRANS = '  CblasNoTrans';
                  } else if (TRANS == 'T') {
                     CTRANS = '    CblasTrans';
                  } else {
                     CTRANS = 'CblasConjTrans';
                  }

                  for (ICD = 1; ICD <= 2; ICD++) { // 70
                     DIAG = ICHD( ICD: ICD );
                     if (DIAG == 'N') {
                        CDIAG = '  CblasNonUnit';
                     } else {
                        CDIAG = '     CblasUnit';
                     }

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 8: 9 ), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

                     for (IX = 1; IX <= NINC; IX++) { // 60
                        INCX = INC( IX );
                        LX = ( INCX ).abs()*N;

                        // Generate the vector X.

                        TRANSL = HALF;
                        zmake('ge', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
                        if ( N > 1 ) {
                           X( N/2 ) = ZERO;
                           XX( 1 + ( INCX ).abs()*( N/2 - 1 ) ) = ZERO;
                        }

                        NC = NC + 1;

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO;
                        TRANSS = TRANS;
                        DIAGS = DIAG;
                        NS = N;
                        KS = K;
                        for (I = 1; I <= LAA; I++) { // 20
                           AS( I ) = AA( I );
                        } // 20
                        LDAS = LDA;
                        for (I = 1; I <= LX; I++) { // 30
                           XS( I ) = XX( I );
                        } // 30
                        INCXS = INCX;

                        // Call the subroutine.

                        if ( SNAME( 10: 11 ) == 'mv' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              cztrmv(IORDER, UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              cztbmv(IORDER, UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              cztpmv(IORDER, UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        } else if ( SNAME( 10: 11 ) == 'sv' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              cztrsv(IORDER, UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              cztbsv(IORDER, UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              cztpsv(IORDER, UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 120;
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLO == UPLOS;
                        ISAME( 2 ) = TRANS == TRANSS;
                        ISAME( 3 ) = DIAG == DIAGS;
                        ISAME( 4 ) = NS == N;
                        if ( FULL ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA );
                           ISAME( 6 ) = LDAS == LDA;
                           if ( NULL ) {
                              ISAME( 7 ) = LZE( XS, XX, LX );
                           } else {
                              ISAME( 7 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME( 8 ) = INCXS == INCX;
                        } else if ( BANDED ) {
                           ISAME( 5 ) = KS == K;
                           ISAME( 6 ) = LZE( AS, AA, LAA );
                           ISAME( 7 ) = LDAS == LDA;
                           if ( NULL ) {
                              ISAME( 8 ) = LZE( XS, XX, LX );
                           } else {
                              ISAME( 8 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME( 9 ) = INCXS == INCX;
                        } else if ( PACKED ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA );
                           if ( NULL ) {
                              ISAME( 6 ) = LZE( XS, XX, LX );
                           } else {
                              ISAME( 6 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ( INCX ).abs() );
                           }
                           ISAME( 7 ) = INCXS == INCX;
                        }

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                        } // 40
                        if ( !SAME ) {
                           FATAL = true;
                           GO TO 120;
                        }

                        if ( !NULL ) {
                           if ( SNAME( 10: 11 ) == 'mv' ) {

                              // Check the result.

                              zmvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, true );
                           } else if ( SNAME( 10: 11 ) == 'sv' ) {

                              // Compute approximation to original vector.

                              for (I = 1; I <= N; I++) { // 50
                                 Z( I ) = XX( 1 + ( I - 1 )* ( INCX ).abs() )                                  XX( 1 + ( I - 1 )*( INCX ).abs() ) = X( I );
                              } // 50
                              zmvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, false );
                           }
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 120;
                        } else {
                           // Avoid repeating tests with N <= 0.
                           GO TO 110;
                        }

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX;
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX;
      }

      } // 130
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT(1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x, I3, ', AP, ', 'X,', I2, ') .' );
 9994 FORMAT(1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x,  2( I3, ',' ), ' A,', I3, ', X,', I2, ') .' );
 9993 FORMAT( 1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x, I3, ', A,', I3, ', X,', I2, ') .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of ZCHK3.

      }
      void zchk4(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER ) {

// Tests ZGERC and ZGERU.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex         ALPHA, ALS, TRANSL;
      double             ERR, ERRMAX;
      int                I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS, NC, ND, NS;
      bool               CONJ, NULL, RESET, SAME;
      // .. Local Arrays ..
      Complex         W( 1 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZGERC, CZGERU, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Executable Statements ..
      CONJ = SNAME( 11: 11 ) == 'c';
      // Define the number of arguments.
      NARGS = 9;

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN );
         ND = N/2 + 1;

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = max( N - ND, 0 );
            IF( IM == 2 ) M = min( N + ND, NMAX );

            // Set LDA to 1 more than minimum value if room.
            LDA = M;
            if (LDA < NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA > NMAX) GO TO 110;
            LAA = LDA*N;
            NULL = N <= 0 || M <= 0;

            for (IX = 1; IX <= NINC; IX++) { // 100
               INCX = INC( IX );
               LX = ( INCX ).abs()*M;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('ge', ' ', ' ', 1, M, X, 1, XX, ( INCX ).abs(), 0, M - 1, RESET, TRANSL );
               if ( M > 1 ) {
                  X( M/2 ) = ZERO;
                  XX( 1 + ( INCX ).abs()*( M/2 - 1 ) ) = ZERO;
               }

               for (IY = 1; IY <= NINC; IY++) { // 90
                  INCY = INC( IY );
                  LY = ( INCY ).abs()*N;

                  // Generate the vector Y.

                  TRANSL = ZERO;
                  zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     Y( N/2 ) = ZERO;
                     YY( 1 + ( INCY ).abs()*( N/2 - 1 ) ) = ZERO;
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 80
                     ALPHA = ALF( IA );

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 8: 9 ), ' ', ' ', M, N, A, NMAX, AA, LDA, M - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1;

                     // Save every datum before calling the subroutine.

                     MS = M;
                     NS = N;
                     ALS = ALPHA;
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I );
                     } // 10
                     LDAS = LDA;
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I );
                     } // 20
                     INCXS = INCX;
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I );
                     } // 30
                     INCYS = INCY;

                     // Call the subroutine.

                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA;
                     if ( CONJ ) {
                        if (REWI) REWIND NTRA;
                        czgerc(IORDER, M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     } else {
                        if (REWI) REWIND NTRA;
                        czgeru(IORDER, M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( !OK ) {
                        WRITE( NOUT, FMT = 9993 );
                        FATAL = true;
                        GO TO 140;
                     }

                     // See what data changed inside subroutine.

                     ISAME( 1 ) = MS == M;
                     ISAME( 2 ) = NS == N;
                     ISAME( 3 ) = ALS == ALPHA;
                     ISAME( 4 ) = LZE( XS, XX, LX );
                     ISAME( 5 ) = INCXS == INCX;
                     ISAME( 6 ) = LZE( YS, YY, LY );
                     ISAME( 7 ) = INCYS == INCY;
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA );
                     } else {
                        ISAME( 8 ) = LZERES( 'ge', ' ', M, N, AS, AA, LDA );
                     }
                     ISAME( 9 ) = LDAS == LDA;

                    // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME && ISAME( I );
                        if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                     } // 40
                     if ( !SAME ) {
                        FATAL = true;
                        GO TO 140;
                     }

                     if ( !NULL ) {

                        // Check the result column by column.

                        if ( INCX > 0 ) {
                           for (I = 1; I <= M; I++) { // 50
                              Z( I ) = X( I );
                           } // 50
                        } else {
                           for (I = 1; I <= M; I++) { // 60
                              Z( I ) = X( M - I + 1 );
                           } // 60
                        }
                        for (J = 1; J <= N; J++) { // 70
                           if ( INCY > 0 ) {
                              W( 1 ) = Y( J );
                           } else {
                              W( 1 ) = Y( N - J + 1 );
                           }
                           if (CONJ) W( 1 ) = DCONJG( W( 1 ) );
                           zmvch('N', M, 1, ALPHA, Z, NMAX, W, 1, ONE, A( 1, J ), 1, YT, G, AA( 1 + ( J - 1 )*LDA ), EPS, ERR, FATAL, NOUT, true );
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 130;
                        } // 70
                     } else {
                        // Avoid repeating tests with M <= 0 or N <= 0.
                        GO TO 110;
                     }

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 150;

      } // 130
      WRITE( NOUT, FMT = 9995 )J;

      } // 140
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA;

      } // 150
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );
 9994 FORMAT(1X, I6, ': ',A12, '(', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ') .' );
 9993 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of ZCHK4.

      }
      void zchk5(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER ) {

// Tests ZHER and ZHPR.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX );
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex         ALPHA, TRANSL;
      double             ERR, ERRMAX, RALPHA, RALS;
      int               I, IA, IC, IN, INCX, INCXS, IX, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      Complex         W( 1 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHER, CZHPR, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG, MAX, DBLE
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ) == 'e';
      PACKED = SNAME( 9: 9 ) == 'p';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 7;
      } else if ( PACKED ) {
         NARGS = 6;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN );
         // Set LDA to 1 more than minimum value if room.
         LDA = N;
         if (LDA < NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA > NMAX) GO TO 100;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2;
         } else {
            LAA = LDA*N;
         }

         for (IC = 1; IC <= 2; IC++) { // 90
            UPLO = ICH( IC: IC );
            if (UPLO == 'U') {
               CUPLO = '    CblasUpper';
            } else {
               CUPLO = '    CblasLower';
            }
            UPPER = UPLO == 'U';

            for (IX = 1; IX <= NINC; IX++) { // 80
               INCX = INC( IX );
               LX = ( INCX ).abs()*N;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('ge', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
               if ( N > 1 ) {
                  X( N/2 ) = ZERO;
                  XX( 1 + ( INCX ).abs()*( N/2 - 1 ) ) = ZERO;
               }

               for (IA = 1; IA <= NALF; IA++) { // 70
                  RALPHA = DBLE( ALF( IA ) );
                  ALPHA = DCMPLX( RALPHA, RZERO );
                  NULL = N <= 0 || RALPHA == RZERO;

                  // Generate the matrix A.

                  TRANSL = ZERO;
                  zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                  NC = NC + 1;

                  // Save every datum before calling the subroutine.

                  UPLOS = UPLO;
                  NS = N;
                  RALS = RALPHA;
                  for (I = 1; I <= LAA; I++) { // 10
                     AS( I ) = AA( I );
                  } // 10
                  LDAS = LDA;
                  for (I = 1; I <= LX; I++) { // 20
                     XS( I ) = XX( I );
                  } // 20
                  INCXS = INCX;

                  // Call the subroutine.

                  if ( FULL ) {
                     if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, RALPHA, INCX, LDA;
                     if (REWI) REWIND NTRA;
                     czher(IORDER, UPLO, N, RALPHA, XX, INCX, AA, LDA );
                  } else if ( PACKED ) {
                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, RALPHA, INCX;
                     if (REWI) REWIND NTRA;
                     czhpr(IORDER, UPLO, N, RALPHA, XX, INCX, AA );
                  }

                  // Check if error-exit was taken incorrectly.

                  if ( !OK ) {
                     WRITE( NOUT, FMT = 9992 );
                     FATAL = true;
                     GO TO 120;
                  }

                  // See what data changed inside subroutines.

                  ISAME( 1 ) = UPLO == UPLOS;
                  ISAME( 2 ) = NS == N;
                  ISAME( 3 ) = RALS == RALPHA;
                  ISAME( 4 ) = LZE( XS, XX, LX );
                  ISAME( 5 ) = INCXS == INCX;
                  if ( NULL ) {
                     ISAME( 6 ) = LZE( AS, AA, LAA );
                  } else {
                    ISAME( 6 ) = LZERES( SNAME( 8: 9 ), UPLO, N, N, AS, AA, LDA );
                  }
                  if ( !PACKED ) {
                     ISAME( 7 ) = LDAS == LDA;
                  }

                  // If data was incorrectly changed, report and return.

                  SAME = true;
                  for (I = 1; I <= NARGS; I++) { // 30
                     SAME = SAME && ISAME( I );
                     if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                  } // 30
                  if ( !SAME ) {
                     FATAL = true;
                     GO TO 120;
                  }

                  if ( !NULL ) {

                     // Check the result column by column.

                     if ( INCX > 0 ) {
                        for (I = 1; I <= N; I++) { // 40
                           Z( I ) = X( I );
                        } // 40
                     } else {
                        for (I = 1; I <= N; I++) { // 50
                           Z( I ) = X( N - I + 1 );
                        } // 50
                     }
                     JA = 1;
                     for (J = 1; J <= N; J++) { // 60
                        W( 1 ) = DCONJG( Z( J ) );
                        if ( UPPER ) {
                           JJ = 1;
                           LJ = J;
                        } else {
                           JJ = J;
                           LJ = N - J + 1;
                        }
                        zmvch('N', LJ, 1, ALPHA, Z( JJ ), LJ, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
                        if ( FULL ) {
                           if ( UPPER ) {
                              JA = JA + LDA;
                           } else {
                              JA = JA + LDA + 1;
                           }
                        } else {
                           JA = JA + LJ;
                        }
                        ERRMAX = max( ERRMAX, ERR );
                        // If got really bad answer, report and return.
                        if (FATAL) GO TO 110;
                     } // 60
                  } else {
                     // Avoid repeating tests if N <= 0.
                     if (N <= 0) GO TO 100;
                  }

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 110
      WRITE( NOUT, FMT = 9995 )J;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, RALPHA, INCX, LDA;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, RALPHA, INCX;
      }

      } // 130
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );
 9994 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',', F4.1, ', X,', I2, ', AP) .' );
 9993 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',', F4.1, ', X,', I2, ', A,', I3, ') .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CZHK5.

      }
      void zchk6(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER ) {

// Tests ZHER2 and ZHPR2.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, HALF, ONE;
      const              ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX, 2 );
      double                         G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      Complex            ALPHA, ALS, TRANSL;
      double                         ERR, ERRMAX;
      int                I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      Complex         W( 2 );
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHER2, CZHPR2, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      const ICH = 'UL';
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ) == 'e';
      PACKED = SNAME( 9: 9 ) == 'p';
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 9;
      } else if ( PACKED ) {
         NARGS = 8;
      }

      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 140
         N = IDIM( IN );
         // Set LDA to 1 more than minimum value if room.
         LDA = N;
         if (LDA < NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA > NMAX) GO TO 140;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2;
         } else {
            LAA = LDA*N;
         }

         for (IC = 1; IC <= 2; IC++) { // 130
            UPLO = ICH( IC: IC );
            if (UPLO == 'U') {
               CUPLO = '    CblasUpper';
            } else {
               CUPLO = '    CblasLower';
            }
            UPPER = UPLO == 'U';

            for (IX = 1; IX <= NINC; IX++) { // 120
               INCX = INC( IX );
               LX = ( INCX ).abs()*N;

               // Generate the vector X.

               TRANSL = HALF;
               zmake('ge', ' ', ' ', 1, N, X, 1, XX, ( INCX ).abs(), 0, N - 1, RESET, TRANSL );
               if ( N > 1 ) {
                  X( N/2 ) = ZERO;
                  XX( 1 + ( INCX ).abs()*( N/2 - 1 ) ) = ZERO;
               }

               for (IY = 1; IY <= NINC; IY++) { // 110
                  INCY = INC( IY );
                  LY = ( INCY ).abs()*N;

                  // Generate the vector Y.

                  TRANSL = ZERO;
                  zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ( INCY ).abs(), 0, N - 1, RESET, TRANSL );
                  if ( N > 1 ) {
                     Y( N/2 ) = ZERO;
                     YY( 1 + ( INCY ).abs()*( N/2 - 1 ) ) = ZERO;
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 100
                     ALPHA = ALF( IA );
                     NULL = N <= 0 || ALPHA == ZERO;

                     // Generate the matrix A.

                     TRANSL = ZERO;
                     zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1;

                     // Save every datum before calling the subroutine.

                     UPLOS = UPLO;
                     NS = N;
                     ALS = ALPHA;
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I );
                     } // 10
                     LDAS = LDA;
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I );
                     } // 20
                     INCXS = INCX;
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I );
                     } // 30
                     INCYS = INCY;

                     // Call the subroutine.

                     if ( FULL ) {
                        if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY, LDA;
                        if (REWI) REWIND NTRA;
                        czher2(IORDER, UPLO, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     } else if ( PACKED ) {
                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY;
                        if (REWI) REWIND NTRA;
                        czhpr2(IORDER, UPLO, N, ALPHA, XX, INCX, YY, INCY, AA );
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( !OK ) {
                        WRITE( NOUT, FMT = 9992 );
                        FATAL = true;
                        GO TO 160;
                     }

                     // See what data changed inside subroutines.

                     ISAME( 1 ) = UPLO == UPLOS;
                     ISAME( 2 ) = NS == N;
                     ISAME( 3 ) = ALS == ALPHA;
                     ISAME( 4 ) = LZE( XS, XX, LX );
                     ISAME( 5 ) = INCXS == INCX;
                     ISAME( 6 ) = LZE( YS, YY, LY );
                     ISAME( 7 ) = INCYS == INCY;
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA );
                     } else {
                        ISAME( 8 ) = LZERES( SNAME( 8: 9 ), UPLO, N, N, AS, AA, LDA );
                     }
                     if ( !PACKED ) {
                        ISAME( 9 ) = LDAS == LDA;
                     }

                    // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME && ISAME( I );
                        if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                     } // 40
                     if ( !SAME ) {
                        FATAL = true;
                        GO TO 160;
                     }

                     if ( !NULL ) {

                        // Check the result column by column.

                        if ( INCX > 0 ) {
                           for (I = 1; I <= N; I++) { // 50
                              Z( I, 1 ) = X( I );
                           } // 50
                        } else {
                           for (I = 1; I <= N; I++) { // 60
                              Z( I, 1 ) = X( N - I + 1 );
                           } // 60
                        }
                        if ( INCY > 0 ) {
                           for (I = 1; I <= N; I++) { // 70
                              Z( I, 2 ) = Y( I );
                           } // 70
                        } else {
                           for (I = 1; I <= N; I++) { // 80
                              Z( I, 2 ) = Y( N - I + 1 );
                           } // 80
                        }
                        JA = 1;
                        for (J = 1; J <= N; J++) { // 90
                           W( 1 ) = ALPHA*DCONJG( Z( J, 2 ) );
                           W( 2 ) = DCONJG( ALPHA )*DCONJG( Z( J, 1 ) );
                           if ( UPPER ) {
                              JJ = 1;
                              LJ = J;
                           } else {
                              JJ = J;
                              LJ = N - J + 1;
                           }
                           zmvch('N', LJ, 2, ONE, Z( JJ, 1 ), NMAX, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
                           if ( FULL ) {
                              if ( UPPER ) {
                                 JA = JA + LDA;
                              } else {
                                 JA = JA + LDA + 1;
                              }
                           } else {
                              JA = JA + LJ;
                           }
                           ERRMAX = max( ERRMAX, ERR );
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 150;
                        } // 90
                     } else {
                        // Avoid repeating tests with N <= 0.
                        if (N <= 0) GO TO 140;
                     }

                  } // 100

               } // 110

            } // 120

         } // 130

      } // 140

      // Report result.

      if ( ERRMAX < THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 170;

      } // 150
      WRITE( NOUT, FMT = 9995 )J;

      } // 160
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY, LDA;
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY;
      }

      } // 170
      return;

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' );
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' );
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );
 9994 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', AP) .' );
 9993 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ') .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of ZCHK6.

      }
      void zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, MV ) {

// Checks the results of the computational tests.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0, RONE = 1.0 ;
      // .. Scalar Arguments ..
      Complex         ALPHA, BETA;
      double             EPS, ERR;
      int                INCX, INCY, M, N, NMAX, NOUT;
      bool               FATAL, MV;
      String             TRANS;
      // .. Array Arguments ..
      Complex         A( NMAX, * ), X( * ), Y( * ), YT( * ), YY( * );
      double             G( * );
      // .. Local Scalars ..
      Complex         C;
      double             ERRI;
      int                I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
      bool               CTRAN, TRAN;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DIMAG, DCONJG, MAX, DBLE, SQRT
      // .. Statement Functions ..
      double             ABS1;
      // .. Statement Function definitions ..
      ABS1( C ) = ( DBLE( C ) ).abs() + ( DIMAG( C ) ).abs();
      // .. Executable Statements ..
      TRAN = TRANS == 'T';
      CTRAN = TRANS == 'C';
      if ( TRAN || CTRAN ) {
         ML = N;
         NL = M;
      } else {
         ML = M;
         NL = N;
      }
      if ( INCX < 0 ) {
         KX = NL;
         INCXL = -1;
      } else {
         KX = 1;
         INCXL = 1;
      }
      if ( INCY < 0 ) {
         KY = ML;
         INCYL = -1;
      } else {
         KY = 1;
         INCYL = 1;
      }

      // Compute expected result in YT using data in A, X and Y.
      // Compute gauges in G.

      IY = KY;
      for (I = 1; I <= ML; I++) { // 40
         YT( IY ) = ZERO;
         G( IY ) = RZERO;
         JX = KX;
         if ( TRAN ) {
            for (J = 1; J <= NL; J++) { // 10
               YT( IY ) = YT( IY ) + A( J, I )*X( JX );
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 10
         } else if ( CTRAN ) {
            for (J = 1; J <= NL; J++) { // 20
               YT( IY ) = YT( IY ) + DCONJG( A( J, I ) )*X( JX );
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 20
         } else {
            for (J = 1; J <= NL; J++) { // 30
               YT( IY ) = YT( IY ) + A( I, J )*X( JX );
               G( IY ) = G( IY ) + ABS1( A( I, J ) )*ABS1( X( JX ) );
               JX = JX + INCXL;
            } // 30
         }
         YT( IY ) = ALPHA*YT( IY ) + BETA*Y( IY );
         G( IY ) = ABS1( ALPHA )*G( IY ) + ABS1( BETA )*ABS1( Y( IY ) );
         IY = IY + INCYL;
      } // 40

      // Compute the error ratio for this result.

      ERR = ZERO;
      for (I = 1; I <= ML; I++) { // 50
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*( INCY ).abs() ) )/EPS;
         if( G( I ) != RZERO ) ERRI = ERRI/G( I );
         ERR = max( ERR, ERRI );
         if( ERR*sqrt( EPS ) >= RONE ) GO TO 60;
      } // 50
      // If the loop completes, all results are at least half accurate.
      GO TO 80;

      // Report fatal error.

   60 FATAL = true;
      WRITE( NOUT, FMT = 9999 );
      for (I = 1; I <= ML; I++) { // 70
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, YT( I ), YY( 1 + ( I - 1 )*( INCY ).abs() );
         } else {
            WRITE( NOUT, FMT = 9998 )I, YY( 1 + ( I - 1 )*( INCY ).abs() ), YT( I );
         }
      } // 70

      } // 80
      return;

 9999 FORMAT(' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'                       EXPECTED RE', 'SULT                    COMPUTED RESULT' );
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) );

      // End of ZMVCH.

      }
      bool lze(RI, RJ, LR ) {

// Tests if two arrays are identical.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      Complex         RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         if( RI( I ) != RJ( I ) ) GO TO 20;
      } // 10
      LZE = true;
      GO TO 30;
      } // 20
      LZE = false;
   30 return;

      // End of LZE.

      }
      bool lzeres(TYPE, UPLO, M, N, AA, AS, LDA ) {

// Tests if selected elements in two arrays are equal.

// TYPE is 'ge', 'he' or 'hp'.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex         AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO == 'U';
      if ( TYPE == 'ge' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 10
         } // 20
      } else if ( TYPE == 'he' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1;
               IEND = J;
            } else {
               IBEG = J;
               IEND = N;
            }
            for (I = 1; I <= IBEG - 1; I++) { // 30
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 30
            for (I = IEND + 1; I <= LDA; I++) { // 40
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 40
         } // 50
      }

      } // 60
      LZERES = true;
      GO TO 80;
      } // 70
      LZERES = false;
   80 return;

      // End of LZERES.

      }
      Complex FUNCTION ZBEG( RESET );

// Generates complex numbers as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ;
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891;
         MJ = 457;
         I = 7;
         J = 7;
         IC = 0;
         RESET = false;
      }

      // The sequence of values of I or J is bounded between 1 and 999.
      // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I or J = 4 or 8, the period will be 25.
      // If initial I or J = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I or J
      // in 6.

      IC = IC + 1;
   10 I = I*MI;
      J = J*MJ;
      I = I - 1000*( I/1000 );
      J = J - 1000*( J/1000 );
      if ( IC >= 5 ) {
         IC = 0;
         GO TO 10;
      }
      ZBEG = DCMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 );
      return;

      // End of ZBEG.

      }
      double ddiff(X, Y ) {

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.

      // .. Scalar Arguments ..
      double               X, Y;
      // .. Executable Statements ..
      DDIFF = X - Y;
      return;

      // End of DDIFF.

      }
      void zmake(TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL ) {

// Generates values for an M by N matrix A within the bandwidth
// defined by KL and KU.
// Stores the values in the array AA in the data structure required
// by the routine, with unwanted elements set to rogue value.

// TYPE is 'ge', 'gb', 'he', 'hb', 'hp', 'tr', 'tb' OR 'tp'.

// Auxiliary routine for test program for Level 2 Blas.

// -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e10, 1.0e10 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      double             RROGUE;
      const              RROGUE = -1.0e10 ;
      // .. Scalar Arguments ..
      Complex         TRANSL;
      int                KL, KU, LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex         A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      //- Complex         ZBEG;
      // EXTERNAL ZBEG
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, DCONJG, MAX, MIN, DBLE
      // .. Executable Statements ..
      GEN = TYPE( 1: 1 ) == 'g';
      SYM = TYPE( 1: 1 ) == 'h';
      TRI = TYPE( 1: 1 ) == 't';
      UPPER = ( SYM || TRI ) && UPLO == 'U';
      LOWER = ( SYM || TRI ) && UPLO == 'L';
      UNIT = TRI && DIAG == 'U';

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN || ( UPPER && I <= J ) || ( LOWER && I >= J ) ) {
               if( ( I <= J && J - I <= KU ) || ( I >= J && I - J <= KL ) ) {
                  A( I, J ) = ZBEG( RESET ) + TRANSL;
               } else {
                  A( I, J ) = ZERO;
               }
               if ( I != J ) {
                  if ( SYM ) {
                     A( J, I ) = DCONJG( A( I, J ) );
                  } else if ( TRI ) {
                     A( J, I ) = ZERO;
                  }
               }
            }
         } // 10
         if (SYM) A( J, J ) = DCMPLX( DBLE( A( J, J ) ), RZERO );
         if( TRI ) A( J, J ) = A( J, J ) + ONE;
         IF( UNIT ) A( J, J ) = ONE;
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE == 'ge' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA( I + ( J - 1 )*LDA ) = A( I, J );
            } // 30
            for (I = M + 1; I <= LDA; I++) { // 40
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 40
         } // 50
      } else if ( TYPE == 'gb' ) {
         for (J = 1; J <= N; J++) { // 90
            for (I1 = 1; I1 <= KU + 1 - J; I1++) { // 60
               AA( I1 + ( J - 1 )*LDA ) = ROGUE;
            } // 60
            DO 70 I2 = I1, min( KL + KU + 1, KU + 1 + M - J );
               AA( I2 + ( J - 1 )*LDA ) = A( I2 + J - KU - 1, J );
            } // 70
            for (I3 = I2; I3 <= LDA; I3++) { // 80
               AA( I3 + ( J - 1 )*LDA ) = ROGUE;
            } // 80
         } // 90
      } else if ( TYPE == 'he' || TYPE == 'tr' ) {
         for (J = 1; J <= N; J++) { // 130
            if ( UPPER ) {
               IBEG = 1;
               if ( UNIT ) {
                  IEND = J - 1;
               } else {
                  IEND = J;
               }
            } else {
               if ( UNIT ) {
                  IBEG = J + 1;
               } else {
                  IBEG = J;
               }
               IEND = N;
            }
            for (I = 1; I <= IBEG - 1; I++) { // 100
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 100
            for (I = IBEG; I <= IEND; I++) { // 110
               AA( I + ( J - 1 )*LDA ) = A( I, J );
            } // 110
            for (I = IEND + 1; I <= LDA; I++) { // 120
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 120
            if ( SYM ) {
               JJ = J + ( J - 1 )*LDA;
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE );
            }
         } // 130
      } else if ( TYPE == 'hb' || TYPE == 'tb' ) {
         for (J = 1; J <= N; J++) { // 170
            if ( UPPER ) {
               KK = KL + 1;
               IBEG = max( 1, KL + 2 - J );
               if ( UNIT ) {
                  IEND = KL;
               } else {
                  IEND = KL + 1;
               }
            } else {
               KK = 1;
               if ( UNIT ) {
                  IBEG = 2;
               } else {
                  IBEG = 1;
               }
               IEND = min( KL + 1, 1 + M - J );
            }
            for (I = 1; I <= IBEG - 1; I++) { // 140
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 140
            for (I = IBEG; I <= IEND; I++) { // 150
               AA( I + ( J - 1 )*LDA ) = A( I + J - KK, J );
            } // 150
            for (I = IEND + 1; I <= LDA; I++) { // 160
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 160
            if ( SYM ) {
               JJ = KK + ( J - 1 )*LDA;
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE );
            }
         } // 170
      } else if ( TYPE == 'hp' || TYPE == 'tp' ) {
         IOFF = 0;
         for (J = 1; J <= N; J++) { // 190
            if ( UPPER ) {
               IBEG = 1;
               IEND = J;
            } else {
               IBEG = J;
               IEND = N;
            }
            for (I = IBEG; I <= IEND; I++) { // 180
               IOFF = IOFF + 1;
               AA( IOFF ) = A( I, J );
               if ( I == J ) {
                  if (UNIT) AA( IOFF ) = ROGUE;
                  IF( SYM ) AA( IOFF ) = DCMPLX( DBLE( AA( IOFF ) ), RROGUE );
               }
            } // 180
         } // 190
      }
      return;

      // End of ZMAKE.

      }
