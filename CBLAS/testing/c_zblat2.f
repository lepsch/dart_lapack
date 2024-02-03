      PROGRAM ZBLAT2

*  Test program for the COMPLEX*16          Level 2 Blas.

*  The program must be driven by a short data file. The first 17 records
*  of the file are read using list-directed input, the last 17 records
*  are read using the format ( A12, L2 ). An annotated example of a data
*  file can be obtained by deleting the first 3 characters from the
*  following 34 lines:
*  'CBLAT2.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
*  -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
*  F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
*  F        LOGICAL FLAG, T TO STOP ON FAILURES.
*  T        LOGICAL FLAG, T TO TEST ERROR EXITS.
*  2        0 TO TEST COLUMN-MAJOR, 1 TO TEST ROW-MAJOR, 2 TO TEST BOTH
*  16.0     THRESHOLD VALUE OF TEST RATIO
*  6                 NUMBER OF VALUES OF N
*  0 1 2 3 5 9       VALUES OF N
*  4                 NUMBER OF VALUES OF K
*  0 1 2 4           VALUES OF K
*  4                 NUMBER OF VALUES OF INCX AND INCY
*  1 2 -1 -2         VALUES OF INCX AND INCY
*  3                 NUMBER OF VALUES OF ALPHA
*  (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
*  3                 NUMBER OF VALUES OF BETA
*  (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
*  cblas_zgemv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zgbmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zhemv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zhbmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zhpmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztrmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztbmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztpmv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztrsv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztbsv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ztpsv  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zgerc  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zgeru  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zher   T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zhpr   T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zher2  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_zhpr2  T PUT F FOR NO TEST. SAME COLUMNS.

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


*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                NSUBS;
      const              NSUBS = 17 ;
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO, RHALF, RONE;
      const              RZERO = 0.0D0, RHALF = 0.5D0, RONE = 1.0D0 ;
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
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BET( NBEMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( 2*NMAX )
      double             G( NMAX );
      int                IDIM( NIDMAX ), INC( NINMAX ), KB( NKBMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      double             DDIFF;
      bool               LZE;
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
      COMMON             /INFOC/INFOT, NOUTC, OK
      COMMON             /SRNAMC/SRNAMT
      // .. Data statements ..
      DATA               SNAMES/'cblas_zgemv ', 'cblas_zgbmv ', 'cblas_zhemv ','cblas_zhbmv ','cblas_zhpmv ', 'cblas_ztrmv ','cblas_ztbmv ','cblas_ztpmv ', 'cblas_ztrsv ','cblas_ztbsv ','cblas_ztpsv ', 'cblas_zgerc ','cblas_zgeru ','cblas_zher  ', 'cblas_zhpr  ','cblas_zher2 ','cblas_zhpr2 '/
      // .. Executable Statements ..

      NOUTC = NOUT

      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SNAPS
      READ( NIN, FMT = * )NTRA
      TRACE = NTRA.GE.0
      if ( TRACE ) {
         OPEN( NTRA, FILE = SNAPS )
      }
      // Read the flag that directs rewinding of the snapshot file.
      READ( NIN, FMT = * )REWI
      REWI = REWI.AND.TRACE
      // Read the flag that directs stopping on any failure.
      READ( NIN, FMT = * )SFATAL
      // Read the flag that indicates whether error exits are to be tested.
      READ( NIN, FMT = * )TSTERR
      // Read the flag that indicates whether row-major data layout to be tested.
      READ( NIN, FMT = * )LAYOUT
      // Read the threshold value of the test ratio
      READ( NIN, FMT = * )THRESH

      // Read and check the parameter values for the tests.

      // Values of N
      READ( NIN, FMT = * )NIDIM
      if ( NIDIM.LT.1.OR.NIDIM.GT.NIDMAX ) {
         WRITE( NOUT, FMT = 9997 )'N', NIDMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM )
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ).LT.0.OR.IDIM( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX
            GO TO 230
         }
      } // 10
      // Values of K
      READ( NIN, FMT = * )NKB
      if ( NKB.LT.1.OR.NKB.GT.NKBMAX ) {
         WRITE( NOUT, FMT = 9997 )'K', NKBMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( KB( I ), I = 1, NKB )
      for (I = 1; I <= NKB; I++) { // 20
         if ( KB( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9995 )
            GO TO 230
         }
      } // 20
      // Values of INCX and INCY
      READ( NIN, FMT = * )NINC
      if ( NINC.LT.1.OR.NINC.GT.NINMAX ) {
         WRITE( NOUT, FMT = 9997 )'INCX AND INCY', NINMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( INC( I ), I = 1, NINC )
      for (I = 1; I <= NINC; I++) { // 30
         if ( INC( I ).EQ.0.OR.ABS( INC( I ) ).GT.INCMAX ) {
            WRITE( NOUT, FMT = 9994 )INCMAX
            GO TO 230
         }
      } // 30
      // Values of ALPHA
      READ( NIN, FMT = * )NALF
      if ( NALF.LT.1.OR.NALF.GT.NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF )
      // Values of BETA
      READ( NIN, FMT = * )NBET
      if ( NBET.LT.1.OR.NBET.GT.NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET )

      // Report values of parameters.

      WRITE( NOUT, FMT = 9993 )
      WRITE( NOUT, FMT = 9992 )( IDIM( I ), I = 1, NIDIM )
      WRITE( NOUT, FMT = 9991 )( KB( I ), I = 1, NKB )
      WRITE( NOUT, FMT = 9990 )( INC( I ), I = 1, NINC )
      WRITE( NOUT, FMT = 9989 )( ALF( I ), I = 1, NALF )
      WRITE( NOUT, FMT = 9988 )( BET( I ), I = 1, NBET )
      if ( .NOT.TSTERR ) {
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9980 )
      }
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9999 )THRESH
      WRITE( NOUT, FMT = * )
      RORDER = .FALSE.
      CORDER = .FALSE.
      if (LAYOUT.EQ.2) {
         RORDER = .TRUE.
         CORDER = .TRUE.
         WRITE( *, FMT = 10002 )
      } else if (LAYOUT.EQ.1) {
         RORDER = .TRUE.
         WRITE( *, FMT = 10001 )
      } else if (LAYOUT.EQ.0) {
         CORDER = .TRUE.
         WRITE( *, FMT = 10000 )
      }
      WRITE( *, FMT = * )

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      for (I = 1; I <= NSUBS; I++) { // 40
         LTEST( I ) = .FALSE.
      } // 40
   50 READ( NIN, FMT = 9984, END = 80 )SNAMET, LTESTT
      for (I = 1; I <= NSUBS; I++) { // 60
         IF( SNAMET.EQ.SNAMES( I ) ) GO TO 70
      } // 60
      WRITE( NOUT, FMT = 9986 )SNAMET
      STOP
   70 LTEST( I ) = LTESTT
      GO TO 50

      } // 80
      CLOSE ( NIN )

      // Compute EPS (the machine precision).

      EPS = RONE
      } // 90
      IF( DDIFF( RONE + EPS, RONE ).EQ.RZERO ) GO TO 100
      EPS = RHALF*EPS
      GO TO 90
      } // 100
      EPS = EPS + EPS
      WRITE( NOUT, FMT = 9998 )EPS

      // Check the reliability of ZMVCH using exact data.

      N = MIN( 32, NMAX )
      for (J = 1; J <= N; J++) { // 120
         for (I = 1; I <= N; I++) { // 110
            A( I, J ) = MAX( I - J + 1, 0 )
         } // 110
         X( J ) = J
         Y( J ) = ZERO
      } // 120
      for (J = 1; J <= N; J++) { // 130
         YY( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
      } // 130
      // YY holds the exact result. On exit from CMVCH YT holds
      // the result computed by CMVCH.
      TRANS = 'N'
      zmvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LZE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }
      TRANS = 'T'
      zmvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LZE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 210
         WRITE( NOUT, FMT = * )
         if ( .NOT.LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM )
         } else {
            SRNAMT = SNAMES( ISNUM )
            // Test error exits.
            if ( TSTERR ) {
               cz2chke(SNAMES( ISNUM ) );
               WRITE( NOUT, FMT = * )
            }
            // Test computations.
            INFOT = 0
            OK = .TRUE.
            FATAL = .FALSE.
            GO TO ( 140, 140, 150, 150, 150, 160, 160, 160, 160, 160, 160, 170, 170, 180, 180, 190, 190 )ISNUM
            // Test ZGEMV, 01, and ZGBMV, 02.
  140       IF (CORDER) THEN
            zchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 0 );
            }
            if (RORDER) {
            zchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 1 );
            }
            GO TO 200
            // Test ZHEMV, 03, ZHBMV, 04, and ZHPMV, 05.
  150      IF (CORDER) THEN
           zchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 0 );
           }
           if (RORDER) {
           zchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, 1 );
           }
            GO TO 200
            // Test ZTRMV, 06, ZTBMV, 07, ZTPMV, 08,
            // ZTRSV, 09, ZTBSV, 10, and ZTPSV, 11.
  160      IF (CORDER) THEN
           zchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200
            // Test ZGERC, 12, ZGERU, 13.
  170      IF (CORDER) THEN
           zchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200
            // Test ZHER, 14, and ZHPR, 15.
  180      IF (CORDER) THEN
           zchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }
            GO TO 200
            // Test ZHER2, 16, and ZHPR2, 17.
  190      IF (CORDER) THEN
           zchk6(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 0 );
           }
           if (RORDER) {
           zchk6(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, 1 );
           }

  200       IF( FATAL.AND.SFATAL ) GO TO 220
         }
      } // 210
      WRITE( NOUT, FMT = 9982 )
      GO TO 240

      } // 220
      WRITE( NOUT, FMT = 9981 )
      GO TO 240

      } // 230
      WRITE( NOUT, FMT = 9987 )

      } // 240
      IF( TRACE ) CLOSE ( NTRA )
      CLOSE ( NOUT )
      STOP

10002 FORMAT( ' COLUMN-MAJOR AND ROW-MAJOR DATA LAYOUTS ARE TESTED' )
10001 FORMAT( ' ROW-MAJOR DATA LAYOUT IS TESTED' )
10000 FORMAT( ' COLUMN-MAJOR DATA LAYOUT IS TESTED' )
 9999 FORMAT(' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 )
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 )
 9997 FORMAT(' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 )
 9995 FORMAT( ' VALUE OF K IS LESS THAN 0' )
 9994 FORMAT( ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ', I2 )
 9993 FORMAT(' TESTS OF THE COMPLEX*16      LEVEL 2 BLAS', //' THE F', 'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9992 FORMAT( '   FOR N              ', 9I6 )
 9991 FORMAT( '   FOR K              ', 7I6 )
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 )
 9989 FORMAT( '   FOR ALPHA          ', 7('(', F4.1, ',', F4.1, ')  ', : ) )
 9988 FORMAT( '   FOR BETA           ', 7('(', F4.1, ',', F4.1, ')  ', : ) )
 9987 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' )
 9986 FORMAT(' SUBPROGRAM NAME ',A12, ' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' )
 9985 FORMAT(' ERROR IN CMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' CMVCH WAS CALLED WITH TRANS = ', A1, ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '.', / ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.' , /' ******* TESTS ABANDONED *******' )
 9984 FORMAT( A12, L2 )
 9983 FORMAT( 1X,A12, ' WAS NOT TESTED' )
 9982 FORMAT( /' END OF TESTS' )
 9981 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' )
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )

      // End of ZBLAT2.

      }
      SUBROUTINE ZCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, IORDER )

*  Tests CGEMV and CGBMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16        ZERO, HALF
      const             ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ) ;
      double            RZERO;
      const             RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, BETA, BLS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA, LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK, NL, NS;
      bool               BANDED, FULL, NULL, RESET, SAME, TRAN;
      String             TRANS, TRANSS;
      String             CTRANS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZGBMV, CZGEMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      DATA               ICH/'NTC'/
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ).EQ.'e'
      BANDED = SNAME( 9: 9 ).EQ.'b'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 11
      } else if ( BANDED ) {
         NARGS = 13
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN )
         ND = N/2 + 1

         for (IM = 1; IM <= 2; IM++) { // 110
            IF( IM.EQ.1 ) M = MAX( N - ND, 0 )             IF( IM.EQ.2 ) M = MIN( N + ND, NMAX )

            if ( BANDED ) {
               NK = NKB
            } else {
               NK = 1
            }
            for (IKU = 1; IKU <= NK; IKU++) { // 100
               if ( BANDED ) {
                  KU = KB( IKU )
                  KL = MAX( KU - 1, 0 )
               } else {
                  KU = N - 1
                  KL = M - 1
               }
               // Set LDA to 1 more than minimum value if room.
               if ( BANDED ) {
                  LDA = KL + KU + 1
               } else {
                  LDA = M
               }
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 100
               LAA = LDA*N
               NULL = N.LE.0.OR.M.LE.0

               // Generate the matrix A.

               TRANSL = ZERO
               zmake(SNAME( 8: 9 ), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL );

               for (IC = 1; IC <= 3; IC++) { // 90
                  TRANS = ICH( IC: IC )
                  if (TRANS.EQ.'N') {
                     CTRANS = '  CblasNoTrans'
                  } else if (TRANS.EQ.'T') {
                     CTRANS = '    CblasTrans'
                  } else {
                     CTRANS = 'CblasConjTrans'
                  }
                  TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'

                  if ( TRAN ) {
                     ML = N
                     NL = M
                  } else {
                     ML = M
                     NL = N
                  }

                  for (IX = 1; IX <= NINC; IX++) { // 80
                     INCX = INC( IX )
                     LX = ABS( INCX )*NL

                     // Generate the vector X.

                     TRANSL = HALF
                     zmake('ge', ' ', ' ', 1, NL, X, 1, XX, ABS( INCX ), 0, NL - 1, RESET, TRANSL );
                     if ( NL.GT.1 ) {
                        X( NL/2 ) = ZERO
                        XX( 1 + ABS( INCX )*( NL/2 - 1 ) ) = ZERO
                     }

                     for (IY = 1; IY <= NINC; IY++) { // 70
                        INCY = INC( IY )
                        LY = ABS( INCY )*ML

                        for (IA = 1; IA <= NALF; IA++) { // 60
                           ALPHA = ALF( IA )

                           for (IB = 1; IB <= NBET; IB++) { // 50
                              BETA = BET( IB )

                              // Generate the vector Y.

                              TRANSL = ZERO
                              zmake('ge', ' ', ' ', 1, ML, Y, 1, YY, ABS( INCY ), 0, ML - 1, RESET, TRANSL );

                              NC = NC + 1

                              // Save every datum before calling the
                              // subroutine.

                              TRANSS = TRANS
                              MS = M
                              NS = N
                              KLS = KL
                              KUS = KU
                              ALS = ALPHA
                              for (I = 1; I <= LAA; I++) { // 10
                                 AS( I ) = AA( I )
                              } // 10
                              LDAS = LDA
                              for (I = 1; I <= LX; I++) { // 20
                                 XS( I ) = XX( I )
                              } // 20
                              INCXS = INCX
                              BLS = BETA
                              for (I = 1; I <= LY; I++) { // 30
                                 YS( I ) = YY( I )
                              } // 30
                              INCYS = INCY

                              // Call the subroutine.

                              if ( FULL ) {
                                 IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CTRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
                                 IF( REWI ) REWIND NTRA                                  CALL CZGEMV( IORDER, TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                              } else if ( BANDED ) {
                                 IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, CTRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
                                 IF( REWI ) REWIND NTRA                                  CALL CZGBMV( IORDER, TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                              }

                             // Check if error-exit was taken incorrectly.

                              if ( .NOT.OK ) {
                                 WRITE( NOUT, FMT = 9993 )
                                 FATAL = .TRUE.
                                 GO TO 130
                              }

                              // See what data changed inside subroutines.

         // IF(TRANS .NE. 'C' .OR. (INCX .GT. 0 .AND. INCY .GT. 0)) THEN
                              ISAME( 1 ) = TRANS.EQ.TRANSS
                              ISAME( 2 ) = MS.EQ.M
                              ISAME( 3 ) = NS.EQ.N
                              if ( FULL ) {
                                 ISAME( 4 ) = ALS.EQ.ALPHA
                                 ISAME( 5 ) = LZE( AS, AA, LAA )
                                 ISAME( 6 ) = LDAS.EQ.LDA
                                 ISAME( 7 ) = LZE( XS, XX, LX )
                                 ISAME( 8 ) = INCXS.EQ.INCX
                                 ISAME( 9 ) = BLS.EQ.BETA
                                 if ( NULL ) {
                                    ISAME( 10 ) = LZE( YS, YY, LY )
                                 } else {
                                    ISAME( 10 ) = LZERES( 'ge', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 11 ) = INCYS.EQ.INCY
                              } else if ( BANDED ) {
                                 ISAME( 4 ) = KLS.EQ.KL
                                 ISAME( 5 ) = KUS.EQ.KU
                                 ISAME( 6 ) = ALS.EQ.ALPHA
                                 ISAME( 7 ) = LZE( AS, AA, LAA )
                                 ISAME( 8 ) = LDAS.EQ.LDA
                                 ISAME( 9 ) = LZE( XS, XX, LX )
                                 ISAME( 10 ) = INCXS.EQ.INCX
                                 ISAME( 11 ) = BLS.EQ.BETA
                                 if ( NULL ) {
                                    ISAME( 12 ) = LZE( YS, YY, LY )
                                 } else {
                                    ISAME( 12 ) = LZERES( 'ge', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 13 ) = INCYS.EQ.INCY
                              }

                              // If data was incorrectly changed, report
                              // and return.

                              SAME = .TRUE.
                              for (I = 1; I <= NARGS; I++) { // 40
                                 SAME = SAME.AND.ISAME( I )
                                 IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                              } // 40
                              if ( .NOT.SAME ) {
                                 FATAL = .TRUE.
                                 GO TO 130
                              }

                              if ( .NOT.NULL ) {

                                 // Check the result.

                                 zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
                                 ERRMAX = MAX( ERRMAX, ERR )
                                 // If got really bad answer, report and
                                 // return.
                                 IF( FATAL ) GO TO 130
                              } else {
                                 // Avoid repeating tests with M.le.0 or
                                 // N.le.0.
                                 GO TO 110
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 140

      } // 130
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CTRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CTRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
      }

      } // 140
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 4( I3, ',' ), '(', F4.1, ',', F4.1, '), A,',/ 10x, I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' )
 9994 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,',/ 10x, I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' )
 9993 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK1.

      }
      SUBROUTINE ZCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, IORDER )

*  Tests CHEMV, CHBMV and CHPMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, BETA, BLS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IK, IN, INCX, INCXS, INCY, INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHBMV, CZHEMV, CZHPMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ).EQ.'e'
      BANDED = SNAME( 9: 9 ).EQ.'b'
      PACKED = SNAME( 9: 9 ).EQ.'p'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 10
      } else if ( BANDED ) {
         NARGS = 11
      } else if ( PACKED ) {
         NARGS = 9
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK )
            } else {
               K = N - 1
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1
            } else {
               LDA = N
            }
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 100
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            for (IC = 1; IC <= 2; IC++) { // 90
               UPLO = ICH( IC: IC )
               if (UPLO.EQ.'U') {
                  CUPLO = '    CblasUpper'
               } else {
                  CUPLO = '    CblasLower'
               }

               // Generate the matrix A.

               TRANSL = ZERO
               zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

               for (IX = 1; IX <= NINC; IX++) { // 80
                  INCX = INC( IX )
                  LX = ABS( INCX )*N

                  // Generate the vector X.

                  TRANSL = HALF
                  zmake('ge', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     X( N/2 ) = ZERO
                     XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IY = 1; IY <= NINC; IY++) { // 70
                     INCY = INC( IY )
                     LY = ABS( INCY )*N

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA )

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB )

                           // Generate the vector Y.

                           TRANSL = ZERO
                           zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           UPLOS = UPLO
                           NS = N
                           KS = K
                           ALS = ALPHA
                           for (I = 1; I <= LAA; I++) { // 10
                              AS( I ) = AA( I )
                           } // 10
                           LDAS = LDA
                           for (I = 1; I <= LX; I++) { // 20
                              XS( I ) = XX( I )
                           } // 20
                           INCXS = INCX
                           BLS = BETA
                           for (I = 1; I <= LY; I++) { // 30
                              YS( I ) = YY( I )
                           } // 30
                           INCYS = INCY

                           // Call the subroutine.

                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, LDA, INCX, BETA, INCY                               IF( REWI ) REWIND NTRA                               CALL CZHEMV( IORDER, UPLO, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, K, ALPHA, LDA, INCX, BETA, INCY
                              IF( REWI ) REWIND NTRA                               CALL CZHBMV( IORDER, UPLO, N, K, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, N, ALPHA, INCX, BETA, INCY
                              IF( REWI ) REWIND NTRA                               CALL CZHPMV( IORDER, UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY )
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( .NOT.OK ) {
                              WRITE( NOUT, FMT = 9992 )
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = UPLO.EQ.UPLOS
                           ISAME( 2 ) = NS.EQ.N
                           if ( FULL ) {
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LZE( AS, AA, LAA )
                              ISAME( 5 ) = LDAS.EQ.LDA
                              ISAME( 6 ) = LZE( XS, XX, LX )
                              ISAME( 7 ) = INCXS.EQ.INCX
                              ISAME( 8 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 9 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 9 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 10 ) = INCYS.EQ.INCY
                           } else if ( BANDED ) {
                              ISAME( 3 ) = KS.EQ.K
                              ISAME( 4 ) = ALS.EQ.ALPHA
                              ISAME( 5 ) = LZE( AS, AA, LAA )
                              ISAME( 6 ) = LDAS.EQ.LDA
                              ISAME( 7 ) = LZE( XS, XX, LX )
                              ISAME( 8 ) = INCXS.EQ.INCX
                              ISAME( 9 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 10 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 10 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 11 ) = INCYS.EQ.INCY
                           } else if ( PACKED ) {
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LZE( AS, AA, LAA )
                              ISAME( 5 ) = LZE( XS, XX, LX )
                              ISAME( 6 ) = INCXS.EQ.INCX
                              ISAME( 7 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 8 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 8 ) = LZERES( 'ge', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 9 ) = INCYS.EQ.INCY
                           }

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = .TRUE.
                           for (I = 1; I <= NARGS; I++) { // 40
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                           } // 40
                           if ( .NOT.SAME ) {
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           if ( .NOT.NULL ) {

                              // Check the result.

                              zmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 120
                           } else {
                              // Avoid repeating tests with N.le.0
                              GO TO 110
                           }

                        } // 50

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, K, ALPHA, LDA, INCX, BETA, INCY
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CUPLO, N, ALPHA, INCX, BETA, INCY
      }

      } // 130
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), AP, X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' )
 9994 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' )
 9993 FORMAT( 1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), A,', I3, ', X,',/ 10x, I2, ',(', F4.1, ',', F4.1, '), ', 'Y,', I2, ') .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of CZHK2.

      }
      SUBROUTINE ZCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, XT, G, Z, IORDER )

*  Tests ZTRMV, ZTBMV, ZTPMV, ZTRSV, ZTBSV and ZTPSV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NIDIM, NINC, NKB, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XT( NMAX ), XX( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         TRANSL
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
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZMAKE, ZMVCH, CZTBMV, CZTBSV, CZTPMV, CZTPSV, CZTRMV, CZTRSV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ).EQ.'r'
      BANDED = SNAME( 9: 9 ).EQ.'b'
      PACKED = SNAME( 9: 9 ).EQ.'p'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 8
      } else if ( BANDED ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 7
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
      // Set up zero vector for ZMVCH.
      for (I = 1; I <= NMAX; I++) { // 10
         Z( I ) = ZERO
      } // 10

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         for (IK = 1; IK <= NK; IK++) { // 100
            if ( BANDED ) {
               K = KB( IK )
            } else {
               K = N - 1
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1
            } else {
               LDA = N
            }
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 100
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            for (ICU = 1; ICU <= 2; ICU++) { // 90
               UPLO = ICHU( ICU: ICU )
               if (UPLO.EQ.'U') {
                  CUPLO = '    CblasUpper'
               } else {
                  CUPLO = '    CblasLower'
               }

               for (ICT = 1; ICT <= 3; ICT++) { // 80
                  TRANS = ICHT( ICT: ICT )
                  if (TRANS.EQ.'N') {
                     CTRANS = '  CblasNoTrans'
                  } else if (TRANS.EQ.'T') {
                     CTRANS = '    CblasTrans'
                  } else {
                     CTRANS = 'CblasConjTrans'
                  }

                  for (ICD = 1; ICD <= 2; ICD++) { // 70
                     DIAG = ICHD( ICD: ICD )
                     if (DIAG.EQ.'N') {
                        CDIAG = '  CblasNonUnit'
                     } else {
                        CDIAG = '     CblasUnit'
                     }

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 8: 9 ), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

                     for (IX = 1; IX <= NINC; IX++) { // 60
                        INCX = INC( IX )
                        LX = ABS( INCX )*N

                        // Generate the vector X.

                        TRANSL = HALF
                        zmake('ge', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
                        if ( N.GT.1 ) {
                           X( N/2 ) = ZERO
                           XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                        }

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        DIAGS = DIAG
                        NS = N
                        KS = K
                        for (I = 1; I <= LAA; I++) { // 20
                           AS( I ) = AA( I )
                        } // 20
                        LDAS = LDA
                        for (I = 1; I <= LX; I++) { // 30
                           XS( I ) = XX( I )
                        } // 30
                        INCXS = INCX

                        // Call the subroutine.

                        if ( SNAME( 10: 11 ).EQ.'mv' ) {
                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTRMV( IORDER, UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTBMV( IORDER, UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTPMV( IORDER, UPLO, TRANS, DIAG, N, AA, XX, INCX )
                           }
                        } else if ( SNAME( 10: 11 ).EQ.'sv' ) {
                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTRSV( IORDER, UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTBSV( IORDER, UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX
                              IF( REWI ) REWIND NTRA                               CALL CZTPSV( IORDER, UPLO, TRANS, DIAG, N, AA, XX, INCX )
                           }
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9992 )
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLO.EQ.UPLOS
                        ISAME( 2 ) = TRANS.EQ.TRANSS
                        ISAME( 3 ) = DIAG.EQ.DIAGS
                        ISAME( 4 ) = NS.EQ.N
                        if ( FULL ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA )
                           ISAME( 6 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 7 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 7 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 8 ) = INCXS.EQ.INCX
                        } else if ( BANDED ) {
                           ISAME( 5 ) = KS.EQ.K
                           ISAME( 6 ) = LZE( AS, AA, LAA )
                           ISAME( 7 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 8 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 8 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 9 ) = INCXS.EQ.INCX
                        } else if ( PACKED ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA )
                           if ( NULL ) {
                              ISAME( 6 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 6 ) = LZERES( 'ge', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 7 ) = INCXS.EQ.INCX
                        }

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                        } // 40
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        if ( .NOT.NULL ) {
                           if ( SNAME( 10: 11 ).EQ.'mv' ) {

                              // Check the result.

                              zmvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, .TRUE. );
                           } else if ( SNAME( 10: 11 ).EQ.'sv' ) {

                              // Compute approximation to original vector.

                              for (I = 1; I <= N; I++) { // 50
                                 Z( I ) = XX( 1 + ( I - 1 )* ABS( INCX ) )                                  XX( 1 + ( I - 1 )*ABS( INCX ) ) = X( I )
                              } // 50
                              zmvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, .FALSE. );
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 120
                        } else {
                           // Avoid repeating tests with N.le.0.
                           GO TO 110
                        }

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, LDA, INCX
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, K, LDA, INCX
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, CUPLO, CTRANS, CDIAG, N, INCX
      }

      } // 130
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT(1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x, I3, ', AP, ', 'X,', I2, ') .' )
 9994 FORMAT(1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x,  2( I3, ',' ), ' A,', I3, ', X,', I2, ') .' )
 9993 FORMAT( 1X, I6, ': ',A12, '(', 3( A14, ',' ),/ 10x, I3, ', A,', I3, ', X,', I2, ') .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK3.

      }
      SUBROUTINE ZCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER )

*  Tests ZGERC and ZGERU.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS, NC, ND, NS;
      bool               CONJ, NULL, RESET, SAME;
      // .. Local Arrays ..
      COMPLEX*16         W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZGERC, CZGERU, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Executable Statements ..
      CONJ = SNAME( 11: 11 ).EQ.'c'
      // Define the number of arguments.
      NARGS = 9

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN )
         ND = N/2 + 1

         for (IM = 1; IM <= 2; IM++) { // 110
            IF( IM.EQ.1 ) M = MAX( N - ND, 0 )             IF( IM.EQ.2 ) M = MIN( N + ND, NMAX )

            // Set LDA to 1 more than minimum value if room.
            LDA = M
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 110
            LAA = LDA*N
            NULL = N.LE.0.OR.M.LE.0

            for (IX = 1; IX <= NINC; IX++) { // 100
               INCX = INC( IX )
               LX = ABS( INCX )*M

               // Generate the vector X.

               TRANSL = HALF
               zmake('ge', ' ', ' ', 1, M, X, 1, XX, ABS( INCX ), 0, M - 1, RESET, TRANSL );
               if ( M.GT.1 ) {
                  X( M/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( M/2 - 1 ) ) = ZERO
               }

               for (IY = 1; IY <= NINC; IY++) { // 90
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 80
                     ALPHA = ALF( IA )

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 8: 9 ), ' ', ' ', M, N, A, NMAX, AA, LDA, M - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     MS = M
                     NS = N
                     ALS = ALPHA
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I )
                     } // 10
                     LDAS = LDA
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I )
                     } // 20
                     INCXS = INCX
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I )
                     } // 30
                     INCYS = INCY

                     // Call the subroutine.

                     IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA
                     if ( CONJ ) {
                        IF( REWI ) REWIND NTRA                         CALL CZGERC( IORDER, M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     } else {
                        IF( REWI ) REWIND NTRA                         CALL CZGERU( IORDER, M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9993 )
                        FATAL = .TRUE.
                        GO TO 140
                     }

                     // See what data changed inside subroutine.

                     ISAME( 1 ) = MS.EQ.M
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LZE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LZE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LZERES( 'ge', ' ', M, N, AS, AA, LDA )
                     }
                     ISAME( 9 ) = LDAS.EQ.LDA

                    // If data was incorrectly changed, report and return.

                     SAME = .TRUE.
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                     } // 40
                     if ( .NOT.SAME ) {
                        FATAL = .TRUE.
                        GO TO 140
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           for (I = 1; I <= M; I++) { // 50
                              Z( I ) = X( I )
                           } // 50
                        } else {
                           for (I = 1; I <= M; I++) { // 60
                              Z( I ) = X( M - I + 1 )
                           } // 60
                        }
                        for (J = 1; J <= N; J++) { // 70
                           if ( INCY.GT.0 ) {
                              W( 1 ) = Y( J )
                           } else {
                              W( 1 ) = Y( N - J + 1 )
                           }
                           IF( CONJ ) W( 1 ) = DCONJG( W( 1 ) )                            CALL ZMVCH( 'N', M, 1, ALPHA, Z, NMAX, W, 1, ONE, A( 1, J ), 1, YT, G, AA( 1 + ( J - 1 )*LDA ), EPS, ERR, FATAL, NOUT, .TRUE. )
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 130
                        } // 70
                     } else {
                        // Avoid repeating tests with M.le.0 or N.le.0.
                        GO TO 110
                     }

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 150

      } // 130
      WRITE( NOUT, FMT = 9995 )J

      } // 140
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA

      } // 150
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT(1X, I6, ': ',A12, '(', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ') .' )
 9993 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK4.

      }
      SUBROUTINE ZCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER )

*  Tests ZHER and ZHPR.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, TRANSL
      double             ERR, ERRMAX, RALPHA, RALS;
      int               I, IA, IC, IN, INCX, INCXS, IX, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      COMPLEX*16         W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHER, CZHPR, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG, MAX, DBLE
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ).EQ.'e'
      PACKED = SNAME( 9: 9 ).EQ.'p'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 7
      } else if ( PACKED ) {
         NARGS = 6
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         IF( LDA.LT.NMAX ) LDA = LDA + 1
         // Skip tests if not enough room.
         IF( LDA.GT.NMAX ) GO TO 100
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         for (IC = 1; IC <= 2; IC++) { // 90
            UPLO = ICH( IC: IC )
            if (UPLO.EQ.'U') {
               CUPLO = '    CblasUpper'
            } else {
               CUPLO = '    CblasLower'
            }
            UPPER = UPLO.EQ.'U'

            for (IX = 1; IX <= NINC; IX++) { // 80
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               zmake('ge', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               for (IA = 1; IA <= NALF; IA++) { // 70
                  RALPHA = DBLE( ALF( IA ) )
                  ALPHA = DCMPLX( RALPHA, RZERO )
                  NULL = N.LE.0.OR.RALPHA.EQ.RZERO

                  // Generate the matrix A.

                  TRANSL = ZERO
                  zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                  NC = NC + 1

                  // Save every datum before calling the subroutine.

                  UPLOS = UPLO
                  NS = N
                  RALS = RALPHA
                  for (I = 1; I <= LAA; I++) { // 10
                     AS( I ) = AA( I )
                  } // 10
                  LDAS = LDA
                  for (I = 1; I <= LX; I++) { // 20
                     XS( I ) = XX( I )
                  } // 20
                  INCXS = INCX

                  // Call the subroutine.

                  if ( FULL ) {
                     IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, RALPHA, INCX, LDA
                     IF( REWI ) REWIND NTRA                      CALL CZHER( IORDER, UPLO, N, RALPHA, XX, INCX, AA, LDA )
                  } else if ( PACKED ) {
                     IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, RALPHA, INCX
                     IF( REWI ) REWIND NTRA                      CALL CZHPR( IORDER, UPLO, N, RALPHA, XX, INCX, AA )
                  }

                  // Check if error-exit was taken incorrectly.

                  if ( .NOT.OK ) {
                     WRITE( NOUT, FMT = 9992 )
                     FATAL = .TRUE.
                     GO TO 120
                  }

                  // See what data changed inside subroutines.

                  ISAME( 1 ) = UPLO.EQ.UPLOS
                  ISAME( 2 ) = NS.EQ.N
                  ISAME( 3 ) = RALS.EQ.RALPHA
                  ISAME( 4 ) = LZE( XS, XX, LX )
                  ISAME( 5 ) = INCXS.EQ.INCX
                  if ( NULL ) {
                     ISAME( 6 ) = LZE( AS, AA, LAA )
                  } else {
                    ISAME( 6 ) = LZERES( SNAME( 8: 9 ), UPLO, N, N, AS, AA, LDA )
                  }
                  if ( .NOT.PACKED ) {
                     ISAME( 7 ) = LDAS.EQ.LDA
                  }

                  // If data was incorrectly changed, report and return.

                  SAME = .TRUE.
                  for (I = 1; I <= NARGS; I++) { // 30
                     SAME = SAME.AND.ISAME( I )
                     IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                  } // 30
                  if ( .NOT.SAME ) {
                     FATAL = .TRUE.
                     GO TO 120
                  }

                  if ( .NOT.NULL ) {

                     // Check the result column by column.

                     if ( INCX.GT.0 ) {
                        for (I = 1; I <= N; I++) { // 40
                           Z( I ) = X( I )
                        } // 40
                     } else {
                        for (I = 1; I <= N; I++) { // 50
                           Z( I ) = X( N - I + 1 )
                        } // 50
                     }
                     JA = 1
                     for (J = 1; J <= N; J++) { // 60
                        W( 1 ) = DCONJG( Z( J ) )
                        if ( UPPER ) {
                           JJ = 1
                           LJ = J
                        } else {
                           JJ = J
                           LJ = N - J + 1
                        }
                        zmvch('N', LJ, 1, ALPHA, Z( JJ ), LJ, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, .TRUE. );
                        if ( FULL ) {
                           if ( UPPER ) {
                              JA = JA + LDA
                           } else {
                              JA = JA + LDA + 1
                           }
                        } else {
                           JA = JA + LJ
                        }
                        ERRMAX = MAX( ERRMAX, ERR )
                        // If got really bad answer, report and return.
                        IF( FATAL ) GO TO 110
                     } // 60
                  } else {
                     // Avoid repeating tests if N.le.0.
                     IF( N.LE.0 ) GO TO 100
                  }

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 110
      WRITE( NOUT, FMT = 9995 )J

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, RALPHA, INCX, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, RALPHA, INCX
      }

      } // 130
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',', F4.1, ', X,', I2, ', AP) .' )
 9993 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',', F4.1, ', X,', I2, ', A,', I3, ') .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of CZHK5.

      }
      SUBROUTINE ZCHK6( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z, IORDER )

*  Tests ZHER2 and ZHPR2.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX, 2 )
      double                         G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16            ALPHA, ALS, TRANSL
      double                         ERR, ERRMAX;
      int                I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             CUPLO;
      String             ICH;
      // .. Local Arrays ..
      COMPLEX*16         W( 2 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL CZHER2, CZHPR2, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool                OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 9: 9 ).EQ.'e'
      PACKED = SNAME( 9: 9 ).EQ.'p'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 8
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 140
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         IF( LDA.LT.NMAX ) LDA = LDA + 1
         // Skip tests if not enough room.
         IF( LDA.GT.NMAX ) GO TO 140
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         for (IC = 1; IC <= 2; IC++) { // 130
            UPLO = ICH( IC: IC )
            if (UPLO.EQ.'U') {
               CUPLO = '    CblasUpper'
            } else {
               CUPLO = '    CblasLower'
            }
            UPPER = UPLO.EQ.'U'

            for (IX = 1; IX <= NINC; IX++) { // 120
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               zmake('ge', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               for (IY = 1; IY <= NINC; IY++) { // 110
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  zmake('ge', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 100
                     ALPHA = ALF( IA )
                     NULL = N.LE.0.OR.ALPHA.EQ.ZERO

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 8: 9 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     UPLOS = UPLO
                     NS = N
                     ALS = ALPHA
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I )
                     } // 10
                     LDAS = LDA
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I )
                     } // 20
                     INCXS = INCX
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I )
                     } // 30
                     INCYS = INCY

                     // Call the subroutine.

                     if ( FULL ) {
                        IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY, LDA
                        IF( REWI ) REWIND NTRA                         CALL CZHER2( IORDER, UPLO, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     } else if ( PACKED ) {
                        IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY
                        IF( REWI ) REWIND NTRA                         CALL CZHPR2( IORDER, UPLO, N, ALPHA, XX, INCX, YY, INCY, AA )
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9992 )
                        FATAL = .TRUE.
                        GO TO 160
                     }

                     // See what data changed inside subroutines.

                     ISAME( 1 ) = UPLO.EQ.UPLOS
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LZE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LZE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LZERES( SNAME( 8: 9 ), UPLO, N, N, AS, AA, LDA )
                     }
                     if ( .NOT.PACKED ) {
                        ISAME( 9 ) = LDAS.EQ.LDA
                     }

                    // If data was incorrectly changed, report and return.

                     SAME = .TRUE.
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                     } // 40
                     if ( .NOT.SAME ) {
                        FATAL = .TRUE.
                        GO TO 160
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           for (I = 1; I <= N; I++) { // 50
                              Z( I, 1 ) = X( I )
                           } // 50
                        } else {
                           for (I = 1; I <= N; I++) { // 60
                              Z( I, 1 ) = X( N - I + 1 )
                           } // 60
                        }
                        if ( INCY.GT.0 ) {
                           for (I = 1; I <= N; I++) { // 70
                              Z( I, 2 ) = Y( I )
                           } // 70
                        } else {
                           for (I = 1; I <= N; I++) { // 80
                              Z( I, 2 ) = Y( N - I + 1 )
                           } // 80
                        }
                        JA = 1
                        for (J = 1; J <= N; J++) { // 90
                           W( 1 ) = ALPHA*DCONJG( Z( J, 2 ) )
                           W( 2 ) = DCONJG( ALPHA )*DCONJG( Z( J, 1 ) )
                           if ( UPPER ) {
                              JJ = 1
                              LJ = J
                           } else {
                              JJ = J
                              LJ = N - J + 1
                           }
                           zmvch('N', LJ, 2, ONE, Z( JJ, 1 ), NMAX, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, .TRUE. );
                           if ( FULL ) {
                              if ( UPPER ) {
                                 JA = JA + LDA
                              } else {
                                 JA = JA + LDA + 1
                              }
                           } else {
                              JA = JA + LJ
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 150
                        } // 90
                     } else {
                        // Avoid repeating tests with N.le.0.
                        IF( N.LE.0 ) GO TO 140
                     }

                  } // 100

               } // 110

            } // 120

         } // 130

      } // 140

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 170

      } // 150
      WRITE( NOUT, FMT = 9995 )J

      } // 160
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, CUPLO, N, ALPHA, INCX, INCY
      }

      } // 170
      RETURN

 9999 FORMAT(' ',A12, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT(' ',A12, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ',A12, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', AP) .' )
 9993 FORMAT(1X, I6, ': ',A12, '(', A14, ',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ') .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK6.

      }
      SUBROUTINE ZMVCH( TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, MV )

*  Checks the results of the computational tests.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D0, 0.0D0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0D0, RONE = 1.0D0 ;
      // .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      double             EPS, ERR;
      int                INCX, INCY, M, N, NMAX, NOUT;
      bool               FATAL, MV;
      String             TRANS;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, * ), X( * ), Y( * ), YT( * ), YY( * )
      double             G( * );
      // .. Local Scalars ..
      COMPLEX*16         C
      double             ERRI;
      int                I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
      bool               CTRAN, TRAN;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DIMAG, DCONJG, MAX, DBLE, SQRT
      // .. Statement Functions ..
      double             ABS1;
      // .. Statement Function definitions ..
      ABS1( C ) = ABS( DBLE( C ) ) + ABS( DIMAG( C ) )
      // .. Executable Statements ..
      TRAN = TRANS.EQ.'T'
      CTRAN = TRANS.EQ.'C'
      if ( TRAN.OR.CTRAN ) {
         ML = N
         NL = M
      } else {
         ML = M
         NL = N
      }
      if ( INCX.LT.0 ) {
         KX = NL
         INCXL = -1
      } else {
         KX = 1
         INCXL = 1
      }
      if ( INCY.LT.0 ) {
         KY = ML
         INCYL = -1
      } else {
         KY = 1
         INCYL = 1
      }

      // Compute expected result in YT using data in A, X and Y.
      // Compute gauges in G.

      IY = KY
      for (I = 1; I <= ML; I++) { // 40
         YT( IY ) = ZERO
         G( IY ) = RZERO
         JX = KX
         if ( TRAN ) {
            for (J = 1; J <= NL; J++) { // 10
               YT( IY ) = YT( IY ) + A( J, I )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
            } // 10
         } else if ( CTRAN ) {
            for (J = 1; J <= NL; J++) { // 20
               YT( IY ) = YT( IY ) + DCONJG( A( J, I ) )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
            } // 20
         } else {
            for (J = 1; J <= NL; J++) { // 30
               YT( IY ) = YT( IY ) + A( I, J )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( I, J ) )*ABS1( X( JX ) )
               JX = JX + INCXL
            } // 30
         }
         YT( IY ) = ALPHA*YT( IY ) + BETA*Y( IY )
         G( IY ) = ABS1( ALPHA )*G( IY ) + ABS1( BETA )*ABS1( Y( IY ) )
         IY = IY + INCYL
      } // 40

      // Compute the error ratio for this result.

      ERR = ZERO
      for (I = 1; I <= ML; I++) { // 50
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*ABS( INCY ) ) )/EPS
         IF( G( I ).NE.RZERO ) ERRI = ERRI/G( I )
         ERR = MAX( ERR, ERRI )
         IF( ERR*SQRT( EPS ).GE.RONE ) GO TO 60
      } // 50
      // If the loop completes, all results are at least half accurate.
      GO TO 80

      // Report fatal error.

   60 FATAL = .TRUE.
      WRITE( NOUT, FMT = 9999 )
      for (I = 1; I <= ML; I++) { // 70
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, YT( I ), YY( 1 + ( I - 1 )*ABS( INCY ) )
         } else {
            WRITE( NOUT, FMT = 9998 )I, YY( 1 + ( I - 1 )*ABS( INCY ) ), YT( I )
         }
      } // 70

      } // 80
      RETURN

 9999 FORMAT(' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'                       EXPECTED RE', 'SULT                    COMPUTED RESULT' )
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) )

      // End of ZMVCH.

      }
      bool    FUNCTION LZE( RI, RJ, LR );

*  Tests if two arrays are identical.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      COMPLEX*16         RI( * ), RJ( * )
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         IF( RI( I ).NE.RJ( I ) ) GO TO 20
      } // 10
      LZE = .TRUE.
      GO TO 30
      } // 20
      LZE = .FALSE.
   30 RETURN

      // End of LZE.

      }
      bool    FUNCTION LZERES( TYPE, UPLO, M, N, AA, AS, LDA );

*  Tests if selected elements in two arrays are equal.

*  TYPE is 'ge', 'he' or 'hp'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX*16         AA( LDA, * ), AS( LDA, * )
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO.EQ.'U'
      if ( TYPE.EQ.'ge' ) {
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = M + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 10
         } // 20
      } else if ( TYPE.EQ.'he' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            DO 30 I = 1, IBEG - 1
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 30
            DO 40 I = IEND + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 40
         } // 50
      }

      } // 60
      LZERES = .TRUE.
      GO TO 80
      } // 70
      LZERES = .FALSE.
   80 RETURN

      // End of LZERES.

      }
      COMPLEX*16 FUNCTION ZBEG( RESET )

*  Generates complex numbers as pairs of random numbers uniformly
*  distributed between -0.5 and 0.5.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891
         MJ = 457
         I = 7
         J = 7
         IC = 0
         RESET = .FALSE.
      }

      // The sequence of values of I or J is bounded between 1 and 999.
      // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I or J = 4 or 8, the period will be 25.
      // If initial I or J = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I or J
      // in 6.

      IC = IC + 1
   10 I = I*MI
      J = J*MJ
      I = I - 1000*( I/1000 )
      J = J - 1000*( J/1000 )
      if ( IC.GE.5 ) {
         IC = 0
         GO TO 10
      }
      ZBEG = DCMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 )
      RETURN

      // End of ZBEG.

      }
      double           FUNCTION DDIFF( X, Y );

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.

      // .. Scalar Arguments ..
      double               X, Y;
      // .. Executable Statements ..
      DDIFF = X - Y
      RETURN

      // End of DDIFF.

      }
      SUBROUTINE ZMAKE( TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL )

*  Generates values for an M by N matrix A within the bandwidth
*  defined by KL and KU.
*  Stores the values in the array AA in the data structure required
*  by the routine, with unwanted elements set to rogue value.

*  TYPE is 'ge', 'gb', 'he', 'hb', 'hp', 'tr', 'tb' OR 'tp'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      COMPLEX*16         ROGUE
      const              ROGUE = ( -1.0D10, 1.0D10 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      double             RROGUE;
      const              RROGUE = -1.0D10 ;
      // .. Scalar Arguments ..
      COMPLEX*16         TRANSL
      int                KL, KU, LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, * ), AA( * )
      // .. Local Scalars ..
      int                I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      COMPLEX*16         ZBEG
      // EXTERNAL ZBEG
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, DCONJG, MAX, MIN, DBLE
      // .. Executable Statements ..
      GEN = TYPE( 1: 1 ).EQ.'g'
      SYM = TYPE( 1: 1 ).EQ.'h'
      TRI = TYPE( 1: 1 ).EQ.'t'
      UPPER = ( SYM.OR.TRI ).AND.UPLO.EQ.'U'
      LOWER = ( SYM.OR.TRI ).AND.UPLO.EQ.'L'
      UNIT = TRI.AND.DIAG.EQ.'U'

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) ) THEN                IF( ( I.LE.J.AND.J - I.LE.KU ).OR. ( I.GE.J.AND.I - J.LE.KL ) ) {
                  A( I, J ) = ZBEG( RESET ) + TRANSL
               } else {
                  A( I, J ) = ZERO
               }
               if ( I.NE.J ) {
                  if ( SYM ) {
                     A( J, I ) = DCONJG( A( I, J ) )
                  } else if ( TRI ) {
                     A( J, I ) = ZERO
                  }
               }
            }
         } // 10
         IF( SYM ) A( J, J ) = DCMPLX( DBLE( A( J, J ) ), RZERO )          IF( TRI ) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE.EQ.'ge' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 30
            DO 40 I = M + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 40
         } // 50
      } else if ( TYPE.EQ.'gb' ) {
         for (J = 1; J <= N; J++) { // 90
            DO 60 I1 = 1, KU + 1 - J
               AA( I1 + ( J - 1 )*LDA ) = ROGUE
            } // 60
            DO 70 I2 = I1, MIN( KL + KU + 1, KU + 1 + M - J )
               AA( I2 + ( J - 1 )*LDA ) = A( I2 + J - KU - 1, J )
            } // 70
            for (I3 = I2; I3 <= LDA; I3++) { // 80
               AA( I3 + ( J - 1 )*LDA ) = ROGUE
            } // 80
         } // 90
      } else if ( TYPE.EQ.'he'.OR.TYPE.EQ.'tr' ) {
         for (J = 1; J <= N; J++) { // 130
            if ( UPPER ) {
               IBEG = 1
               if ( UNIT ) {
                  IEND = J - 1
               } else {
                  IEND = J
               }
            } else {
               if ( UNIT ) {
                  IBEG = J + 1
               } else {
                  IBEG = J
               }
               IEND = N
            }
            DO 100 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 100
            for (I = IBEG; I <= IEND; I++) { // 110
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 110
            DO 120 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 120
            if ( SYM ) {
               JJ = J + ( J - 1 )*LDA
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE )
            }
         } // 130
      } else if ( TYPE.EQ.'hb'.OR.TYPE.EQ.'tb' ) {
         for (J = 1; J <= N; J++) { // 170
            if ( UPPER ) {
               KK = KL + 1
               IBEG = MAX( 1, KL + 2 - J )
               if ( UNIT ) {
                  IEND = KL
               } else {
                  IEND = KL + 1
               }
            } else {
               KK = 1
               if ( UNIT ) {
                  IBEG = 2
               } else {
                  IBEG = 1
               }
               IEND = MIN( KL + 1, 1 + M - J )
            }
            DO 140 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 140
            for (I = IBEG; I <= IEND; I++) { // 150
               AA( I + ( J - 1 )*LDA ) = A( I + J - KK, J )
            } // 150
            DO 160 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 160
            if ( SYM ) {
               JJ = KK + ( J - 1 )*LDA
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE )
            }
         } // 170
      } else if ( TYPE.EQ.'hp'.OR.TYPE.EQ.'tp' ) {
         IOFF = 0
         for (J = 1; J <= N; J++) { // 190
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            for (I = IBEG; I <= IEND; I++) { // 180
               IOFF = IOFF + 1
               AA( IOFF ) = A( I, J )
               if ( I.EQ.J ) {
                  IF( UNIT ) AA( IOFF ) = ROGUE                   IF( SYM ) AA( IOFF ) = DCMPLX( DBLE( AA( IOFF ) ), RROGUE )
               }
            } // 180
         } // 190
      }
      RETURN

      // End of ZMAKE.

      }
