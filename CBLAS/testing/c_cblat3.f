      void main() {
// Test program for the COMPLEX          Level 3 Blas.

// The program must be driven by a short data file. The first 13 records
// of the file are read using list-directed input, the last 9 records
// are read using the format ( A12, L2 ). An annotated example of a data
// file can be obtained by deleting the first 3 characters from the
// following 22 lines:
// 'CBLAT3.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
// -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF < 0)
// F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
// F        LOGICAL FLAG, T TO STOP ON FAILURES.
// T        LOGICAL FLAG, T TO TEST ERROR EXITS.
// 2        0 TO TEST COLUMN-MAJOR, 1 TO TEST ROW-MAJOR, 2 TO TEST BOTH
// 16.0     THRESHOLD VALUE OF TEST RATIO
// 6                 NUMBER OF VALUES OF N
// 0 1 2 3 5 9       VALUES OF N
// 3                 NUMBER OF VALUES OF ALPHA
// (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
// 3                 NUMBER OF VALUES OF BETA
// (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
// cblas_cgemm  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_chemm  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_csymm  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ctrmm  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_ctrsm  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_cherk  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_csyrk  T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_cher2k T PUT F FOR NO TEST. SAME COLUMNS.
// cblas_csyr2k T PUT F FOR NO TEST. SAME COLUMNS.

// See:

      // Dongarra J. J., Du Croz J. J., Duff I. S. and Hammarling S.
      // A Set of Level 3 Basic Linear Algebra Subprograms.

      // Technical Memorandum No.88 (Revision 1), Mathematics and
      // Computer Science Division, Argonne National Laboratory, 9700
      // South Cass Avenue, Argonne, Illinois 60439, US.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                NSUBS;
      const              NSUBS = 9 ;
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO, RHALF, RONE;
      const              RZERO = 0.0, RHALF = 0.5, RONE = 1.0 ;
      int                NMAX;
      const              NMAX = 65 ;
      int                NIDMAX, NALMAX, NBEMAX;
      const              NIDMAX = 9, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      REAL               EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NTRA, LAYOUT;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR, CORDER, RORDER;
      String             TRANSA, TRANSB;
      String             SNAMET;
      String             SNAPS;
      // .. Local Arrays ..
      COMPLEX            AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBEMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      REAL               SDIFF;
      bool               LCE;
      // EXTERNAL SDIFF, LCE
      // .. External Subroutines ..
      // EXTERNAL CCHK1, CCHK2, CCHK3, CCHK4, CCHK5, CMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // COMMON /SRNAMC/SRNAMT
      // .. Data statements ..
      DATA               SNAMES/'cblas_cgemm ', 'cblas_chemm ', 'cblas_csymm ', 'cblas_ctrmm ', 'cblas_ctrsm ', 'cblas_cherk ', 'cblas_csyrk ', 'cblas_cher2k', 'cblas_csyr2k'/;
      // .. Executable Statements ..

      NOUTC = NOUT;

      // Read name and unit number for snapshot output file and open file.

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
         GO TO 220;
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM );
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ) < 0 || IDIM( I ) > NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX;
            GO TO 220;
         }
      } // 10
      // Values of ALPHA
      READ( NIN, FMT = * )NALF;
      if ( NALF < 1 || NALF > NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX;
         GO TO 220;
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF );
      // Values of BETA
      READ( NIN, FMT = * )NBET;
      if ( NBET < 1 || NBET > NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX;
         GO TO 220;
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET );

      // Report values of parameters.

      WRITE( NOUT, FMT = 9995 );
      WRITE( NOUT, FMT = 9994 )( IDIM( I ), I = 1, NIDIM );
      WRITE( NOUT, FMT = 9993 )( ALF( I ), I = 1, NALF );
      WRITE( NOUT, FMT = 9992 )( BET( I ), I = 1, NBET );
      if ( !TSTERR ) {
         WRITE( NOUT, FMT = * );
         WRITE( NOUT, FMT = 9984 );
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

      for (I = 1; I <= NSUBS; I++) { // 20
         LTEST( I ) = false;
      } // 20
   30 READ( NIN, FMT = 9988, END = 60 )SNAMET, LTESTT;
      for (I = 1; I <= NSUBS; I++) { // 40
         IF( SNAMET == SNAMES( I ) ) GO TO 50;
      } // 40
      WRITE( NOUT, FMT = 9990 )SNAMET;
      STOP;
   50 LTEST( I ) = LTESTT;
      GO TO 30;

      } // 60
      CLOSE ( NIN );

      // Compute EPS (the machine precision).

      EPS = RONE;
      } // 70
      IF( SDIFF( RONE + EPS, RONE ) == RZERO ) GO TO 80;
      EPS = RHALF*EPS;
      GO TO 70;
      } // 80
      EPS = EPS + EPS;
      WRITE( NOUT, FMT = 9998 )EPS;

      // Check the reliability of CMMCH using exact data.

      N = MIN( 32, NMAX );
      for (J = 1; J <= N; J++) { // 100
         for (I = 1; I <= N; I++) { // 90
            AB( I, J ) = MAX( I - J + 1, 0 );
         } // 90
         AB( J, NMAX + 1 ) = J;
         AB( 1, NMAX + J ) = J;
         C( J, 1 ) = ZERO;
      } // 100
      for (J = 1; J <= N; J++) { // 110
         CC( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
      } // 110
      // CC holds the exact result. On exit from CMMCH CT holds
      // the result computed by CMMCH.
      TRANSA = 'N';
      TRANSB = 'N';
      cmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LCE( CC, CT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }
      TRANSB = 'C';
      cmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LCE( CC, CT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }
      for (J = 1; J <= N; J++) { // 120
         AB( J, NMAX + 1 ) = N - J + 1;
         AB( 1, NMAX + J ) = N - J + 1;
      } // 120
      for (J = 1; J <= N; J++) { // 130
         CC( N - J + 1 ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
      } // 130
      TRANSA = 'C';
      TRANSB = 'N';
      cmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LCE( CC, CT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }
      TRANSB = 'C';
      cmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LCE( CC, CT, N );
      if ( !SAME || ERR != RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 200
         WRITE( NOUT, FMT = * );
         if ( !LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9987 )SNAMES( ISNUM );
         } else {
            SRNAMT = SNAMES( ISNUM );
            // Test error exits.
            if ( TSTERR ) {
               cc3chke(SNAMES( ISNUM ) );
               WRITE( NOUT, FMT = * );
            }
            // Test computations.
            INFOT = 0;
            OK = true;
            FATAL = false;
            GO TO ( 140, 150, 150, 160, 160, 170, 170, 180, 180 )ISNUM;
            // Test CGEMM, 01.
  140       IF (CORDER) THEN;
            cchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 );
            }
            if (RORDER) {
            cchk1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 );
            }
            GO TO 190;
            // Test CHEMM, 02, CSYMM, 03.
  150       IF (CORDER) THEN;
            cchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 );
            }
            if (RORDER) {
            cchk2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 );
            }
            GO TO 190;
            // Test CTRMM, 04, CTRSM, 05.
  160       IF (CORDER) THEN;
            cchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C, 0 );
            }
            if (RORDER) {
            cchk3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C, 1 );
            }
            GO TO 190;
            // Test CHERK, 06, CSYRK, 07.
  170       IF (CORDER) THEN;
            cchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 );
            }
            if (RORDER) {
            cchk4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 );
            }
            GO TO 190;
            // Test CHER2K, 08, CSYR2K, 09.
  180       IF (CORDER) THEN;
            cchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, 0 );
            }
            if (RORDER) {
            cchk5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, 1 );
            }
            GO TO 190;

  190       IF( FATAL && SFATAL ) GO TO 210;
         }
      } // 200
      WRITE( NOUT, FMT = 9986 );
      GO TO 230;

      } // 210
      WRITE( NOUT, FMT = 9985 );
      GO TO 230;

      } // 220
      WRITE( NOUT, FMT = 9991 );

      } // 230
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT );
      STOP;

10002 FORMAT( ' COLUMN-MAJOR AND ROW-MAJOR DATA LAYOUTS ARE TESTED' )
10001 FORMAT(' ROW-MAJOR DATA LAYOUT IS TESTED' )
10000 FORMAT(' COLUMN-MAJOR DATA LAYOUT IS TESTED' )
 9999 FORMAT(' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 );
 9998 FORMAT(' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 );
 9997 FORMAT(' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 );
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 );
 9995 FORMAT(' TESTS OF THE COMPLEX          LEVEL 3 BLAS', //' THE F', 'OLLOWING PARAMETER VALUES WILL BE USED:' );
 9994 FORMAT( '   FOR N              ', 9I6 );
 9993 FORMAT( '   FOR ALPHA          ', 7( '(', F4.1, ',', F4.1, ')  ', : ) );
 9992 FORMAT( '   FOR BETA           ', 7( '(', F4.1, ',', F4.1, ')  ', : ) );
 9991 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' );
 9990 FORMAT(' SUBPROGRAM NAME ', A12,' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' );
 9989 FORMAT(' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' CMMCH WAS CALLED WITH TRANSA = ', A1, 'AND TRANSB = ', A1, /' AND RETURNED SAME = ', L1, ' AND ', ' ERR = ', F12.3, '.', /' THIS MAY BE DUE TO FAULTS IN THE ', 'ARITHMETIC OR THE COMPILER.', /' ******* TESTS ABANDONED ', '*******' );
 9988 FORMAT( A12,L2 );
 9987 FORMAT( 1X, A12,' WAS NOT TESTED' );
 9986 FORMAT( /' END OF TESTS' );
 9985 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' );
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' );

      // End of CBLAT3.

      }
      SUBROUTINE CCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER );

// Tests CGEMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS;
      REAL               ERR, ERRMAX;
      int                I, IA, IB, ICA, ICB, IK, IM, IN, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MA, MB, MS, N, NA, NARGS, NB, NC, NS;
      bool               NULL, RESET, SAME, TRANA, TRANB;
      String             TRANAS, TRANBS, TRANSA, TRANSB;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CCGEMM, CMAKE, CMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'NTC'/;
      // .. Executable Statements ..

      NARGS = 13;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IM = 1; IM <= NIDIM; IM++) { // 110
         M = IDIM( IM );

         for (IN = 1; IN <= NIDIM; IN++) { // 100
            N = IDIM( IN );
            // Set LDC to 1 more than minimum value if room.
            LDC = M;
            if (LDC < NMAX) LDC = LDC + 1;
            // Skip tests if not enough room.
            if (LDC > NMAX) GO TO 100;
            LCC = LDC*N;
            NULL = N <= 0 || M <= 0;

            for (IK = 1; IK <= NIDIM; IK++) { // 90
               K = IDIM( IK );

               for (ICA = 1; ICA <= 3; ICA++) { // 80
                  TRANSA = ICH( ICA: ICA );
                  TRANA = TRANSA == 'T' || TRANSA == 'C';

                  if ( TRANA ) {
                     MA = K;
                     NA = M;
                  } else {
                     MA = M;
                     NA = K;
                  }
                  // Set LDA to 1 more than minimum value if room.
                  LDA = MA;
                  if (LDA < NMAX) LDA = LDA + 1;
                  // Skip tests if not enough room.
                  if (LDA > NMAX) GO TO 80;
                  LAA = LDA*NA;

                  // Generate the matrix A.

                  cmake('ge', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (ICB = 1; ICB <= 3; ICB++) { // 70
                     TRANSB = ICH( ICB: ICB );
                     TRANB = TRANSB == 'T' || TRANSB == 'C';

                     if ( TRANB ) {
                        MB = N;
                        NB = K;
                     } else {
                        MB = K;
                        NB = N;
                     }
                     // Set LDB to 1 more than minimum value if room.
                     LDB = MB;
                     if (LDB < NMAX) LDB = LDB + 1;
                     // Skip tests if not enough room.
                     if (LDB > NMAX) GO TO 70;
                     LBB = LDB*NB;

                     // Generate the matrix B.

                     cmake('ge', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO );

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA );

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB );

                           // Generate the matrix C.

                           cmake('ge', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

                           NC = NC + 1;

                           // Save every datum before calling the
                           // subroutine.

                           TRANAS = TRANSA;
                           TRANBS = TRANSB;
                           MS = M;
                           NS = N;
                           KS = K;
                           ALS = ALPHA;
                           for (I = 1; I <= LAA; I++) { // 10
                              AS( I ) = AA( I );
                           } // 10
                           LDAS = LDA;
                           for (I = 1; I <= LBB; I++) { // 20
                              BS( I ) = BB( I );
                           } // 20
                           LDBS = LDB;
                           BLS = BETA;
                           for (I = 1; I <= LCC; I++) { // 30
                              CS( I ) = CC( I );
                           } // 30
                           LDCS = LDC;

                           // Call the subroutine.

                           if (TRACE) CALL CPRCN1(NTRA, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC);
                           if (REWI) REWIND NTRA;
                           ccgemm(IORDER, TRANSA, TRANSB, M, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                           // Check if error-exit was taken incorrectly.

                           if ( !OK ) {
                              WRITE( NOUT, FMT = 9994 );
                              FATAL = true;
                              GO TO 120;
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = TRANSA == TRANAS;
                           ISAME( 2 ) = TRANSB == TRANBS;
                           ISAME( 3 ) = MS == M;
                           ISAME( 4 ) = NS == N;
                           ISAME( 5 ) = KS == K;
                           ISAME( 6 ) = ALS == ALPHA;
                           ISAME( 7 ) = LCE( AS, AA, LAA );
                           ISAME( 8 ) = LDAS == LDA;
                           ISAME( 9 ) = LCE( BS, BB, LBB );
                           ISAME( 10 ) = LDBS == LDB;
                           ISAME( 11 ) = BLS == BETA;
                           if ( NULL ) {
                              ISAME( 12 ) = LCE( CS, CC, LCC );
                           } else {
                             ISAME( 12 ) = LCERES( 'ge', ' ', M, N, CS, CC, LDC );
                           }
                           ISAME( 13 ) = LDCS == LDC;

                           // If data was incorrectly changed, report
                           // and return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 40
                              SAME = SAME && ISAME( I );
                              IF( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                           } // 40
                           if ( !SAME ) {
                              FATAL = true;
                              GO TO 120;
                           }

                           if ( !NULL ) {

                              // Check the result.

                             cmmch(TRANSA, TRANSB, M, N, K, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
                              ERRMAX = MAX( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 120;
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
         if (IORDER == 0) WRITE( NOUT, FMT = 10000 )SNAME, NC;
         if (IORDER == 1) WRITE( NOUT, FMT = 10001 )SNAME, NC;
      } else {
         if (IORDER == 0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX;
         if (IORDER == 1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      cprcn1(NOUT, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC);

      } // 130
      return;

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' );
 9995 FORMAT( 1X, I6, ': ', A12,'(''', A1, ''',''', A1, ''',', 3( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1, ',', F4.1, '), C,', I3, ').' );
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CCHK1.

      }

      SUBROUTINE CPRCN1(NOUT, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC);
      int              NOUT, NC, IORDER, M, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA;
      String           TRANSA, TRANSB;
      String           SNAME;
      String           CRC, CTA,CTB;

      if (TRANSA == 'N') {
         CTA = '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CTA = '    CblasTrans';
      } else {
         CTA = 'CblasConjTrans';
      }
      if (TRANSB == 'N') {
         CTB = '  CblasNoTrans';
      } else if (TRANSB == 'T') {
         CTB = '    CblasTrans';
      } else {
         CTB = 'CblasConjTrans';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CTA,CTB;
      WRITE(NOUT, FMT = 9994)M, N, K, ALPHA, LDA, LDB, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',');
 9994 FORMAT( 10X, 3( I3, ',' ) ,' (', F4.1,',',F4.1,') , A,', I3, ', B,', I3, ', (', F4.1,',',F4.1,') , C,', I3, ').' );
      }

      SUBROUTINE CCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER );

// Tests CHEMM and CSYMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS;
      REAL               ERR, ERRMAX;
      int                I, IA, IB, ICS, ICU, IM, IN, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MS, N, NA, NARGS, NC, NS;
      bool               CONJ, LEFT, NULL, RESET, SAME;
      String             SIDE, SIDES, UPLO, UPLOS;
      String             ICHS, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CCHEMM, CMAKE, CMMCH, CCSYMM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHS/'LR'/, ICHU/'UL'/;
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ) == 'he';

      NARGS = 12;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IM = 1; IM <= NIDIM; IM++) { // 100
         M = IDIM( IM );

         for (IN = 1; IN <= NIDIM; IN++) { // 90
            N = IDIM( IN );
            // Set LDC to 1 more than minimum value if room.
            LDC = M;
            if (LDC < NMAX) LDC = LDC + 1;
            // Skip tests if not enough room.
            if (LDC > NMAX) GO TO 90;
            LCC = LDC*N;
            NULL = N <= 0 || M <= 0;
            // Set LDB to 1 more than minimum value if room.
            LDB = M;
            if (LDB < NMAX) LDB = LDB + 1;
            // Skip tests if not enough room.
            if (LDB > NMAX) GO TO 90;
            LBB = LDB*N;

            // Generate the matrix B.

            cmake('ge', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

            for (ICS = 1; ICS <= 2; ICS++) { // 80
               SIDE = ICHS( ICS: ICS );
               LEFT = SIDE == 'L';

               if ( LEFT ) {
                  NA = M;
               } else {
                  NA = N;
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA;
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 80;
               LAA = LDA*NA;

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU );

                  // Generate the hermitian or symmetric matrix A.

                  cmake(SNAME( 8: 9 ), UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );

                        // Generate the matrix C.

                        cmake('ge', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1;

                        // Save every datum before calling the
                        // subroutine.

                        SIDES = SIDE;
                        UPLOS = UPLO;
                        MS = M;
                        NS = N;
                        ALS = ALPHA;
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I );
                        } // 10
                        LDAS = LDA;
                        for (I = 1; I <= LBB; I++) { // 20
                           BS( I ) = BB( I );
                        } // 20
                        LDBS = LDB;
                        BLS = BETA;
                        for (I = 1; I <= LCC; I++) { // 30
                           CS( I ) = CC( I );
                        } // 30
                        LDCS = LDC;

                        // Call the subroutine.

                        if (TRACE) CALL CPRCN2(NTRA, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
                        if (REWI) REWIND NTRA;
                        if ( CONJ ) {
                           cchemm(IORDER, SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        } else {
                           ccsymm(IORDER, SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9994 );
                           FATAL = true;
                           GO TO 110;
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = SIDES == SIDE;
                        ISAME( 2 ) = UPLOS == UPLO;
                        ISAME( 3 ) = MS == M;
                        ISAME( 4 ) = NS == N;
                        ISAME( 5 ) = ALS == ALPHA;
                        ISAME( 6 ) = LCE( AS, AA, LAA );
                        ISAME( 7 ) = LDAS == LDA;
                        ISAME( 8 ) = LCE( BS, BB, LBB );
                        ISAME( 9 ) = LDBS == LDB;
                        ISAME( 10 ) = BLS == BETA;
                        if ( NULL ) {
                           ISAME( 11 ) = LCE( CS, CC, LCC );
                        } else {
                           ISAME( 11 ) = LCERES( 'ge', ' ', M, N, CS, CC, LDC );
                        }
                        ISAME( 12 ) = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           IF( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                        } // 40
                        if ( !SAME ) {
                           FATAL = true;
                           GO TO 110;
                        }

                        if ( !NULL ) {

                           // Check the result.

                           if ( LEFT ) {
                              cmmch('N', 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
                           } else {
                              cmmch('N', 'N', M, N, N, ALPHA, B, NMAX, A, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
                           }
                           ERRMAX = MAX( ERRMAX, ERR );
                           // If got really bad answer, report and
                           // return.
                           if (FATAL) GO TO 110;
                        }

                     } // 50

                  } // 60

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX < THRESH ) {
         if (IORDER == 0) WRITE( NOUT, FMT = 10000 )SNAME, NC;
         if (IORDER == 1) WRITE( NOUT, FMT = 10001 )SNAME, NC;
      } else {
         if (IORDER == 0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX;
         if (IORDER == 1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX;
      }
      GO TO 120;

      } // 110
      WRITE( NOUT, FMT = 9996 )SNAME;
      cprcn2(NOUT, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);

      } // 120
      return;

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' );
 9995 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1, ',', F4.1, '), C,', I3, ')    .' );
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CCHK2.

      }

      SUBROUTINE CPRCN2(NOUT, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
      int              NOUT, NC, IORDER, M, N, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA;
      String           SIDE, UPLO;
      String           SNAME;
      String           CRC, CS,CU;

      if (SIDE == 'L') {
         CS =  '     CblasLeft';
      } else {
         CS =  '    CblasRight';
      }
      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CS,CU;
      WRITE(NOUT, FMT = 9994)M, N, ALPHA, LDA, LDB, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',');
 9994 FORMAT( 10X, 2( I3, ',' ),' (',F4.1,',',F4.1, '), A,', I3, ', B,', I3, ', (',F4.1,',',F4.1, '), ', 'C,', I3, ').' );
      }

      SUBROUTINE CCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, A, AA, AS, B, BB, BS, CT, G, C, IORDER );

// Tests CTRMM and CTRSM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH;
      int                NALF, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CT( NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS;
      REAL               ERR, ERRMAX;
      int               I, IA, ICD, ICS, ICT, ICU, IM, IN, J, LAA, LBB, LDA, LDAS, LDB, LDBS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String            DIAG, DIAGS, SIDE, SIDES, TRANAS, TRANSA, UPLO, UPLOS;
      String             ICHD, ICHS, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CMAKE, CMMCH, CCTRMM, CCTRSM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA              ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/, ICHS/'LR'/;
      // .. Executable Statements ..

      NARGS = 11;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;
      // Set up zero matrix for CMMCH.
      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            C( I, J ) = ZERO;
         } // 10
      } // 20

      for (IM = 1; IM <= NIDIM; IM++) { // 140
         M = IDIM( IM );

         for (IN = 1; IN <= NIDIM; IN++) { // 130
            N = IDIM( IN );
            // Set LDB to 1 more than minimum value if room.
            LDB = M;
            if (LDB < NMAX) LDB = LDB + 1;
            // Skip tests if not enough room.
            if (LDB > NMAX) GO TO 130;
            LBB = LDB*N;
            NULL = M <= 0 || N <= 0;

            for (ICS = 1; ICS <= 2; ICS++) { // 120
               SIDE = ICHS( ICS: ICS );
               LEFT = SIDE == 'L';
               if ( LEFT ) {
                  NA = M;
               } else {
                  NA = N;
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA;
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 130;
               LAA = LDA*NA;

               for (ICU = 1; ICU <= 2; ICU++) { // 110
                  UPLO = ICHU( ICU: ICU );

                  for (ICT = 1; ICT <= 3; ICT++) { // 100
                     TRANSA = ICHT( ICT: ICT );

                     for (ICD = 1; ICD <= 2; ICD++) { // 90
                        DIAG = ICHD( ICD: ICD );

                        for (IA = 1; IA <= NALF; IA++) { // 80
                           ALPHA = ALF( IA );

                           // Generate the matrix A.

                           cmake('tr', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                           // Generate the matrix B.

                           cmake('ge', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

                           NC = NC + 1;

                           // Save every datum before calling the
                           // subroutine.

                           SIDES = SIDE;
                           UPLOS = UPLO;
                           TRANAS = TRANSA;
                           DIAGS = DIAG;
                           MS = M;
                           NS = N;
                           ALS = ALPHA;
                           for (I = 1; I <= LAA; I++) { // 30
                              AS( I ) = AA( I );
                           } // 30
                           LDAS = LDA;
                           for (I = 1; I <= LBB; I++) { // 40
                              BS( I ) = BB( I );
                           } // 40
                           LDBS = LDB;

                           // Call the subroutine.

                           if ( SNAME( 10: 11 ) == 'mm' ) {
                              if (TRACE) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);
                              if (REWI) REWIND NTRA;
                              cctrmm(IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           } else if ( SNAME( 10: 11 ) == 'sm' ) {
                              if (TRACE) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);
                              if (REWI) REWIND NTRA;
                              cctrsm(IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( !OK ) {
                              WRITE( NOUT, FMT = 9994 );
                              FATAL = true;
                              GO TO 150;
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = SIDES == SIDE;
                           ISAME( 2 ) = UPLOS == UPLO;
                           ISAME( 3 ) = TRANAS == TRANSA;
                           ISAME( 4 ) = DIAGS == DIAG;
                           ISAME( 5 ) = MS == M;
                           ISAME( 6 ) = NS == N;
                           ISAME( 7 ) = ALS == ALPHA;
                           ISAME( 8 ) = LCE( AS, AA, LAA );
                           ISAME( 9 ) = LDAS == LDA;
                           if ( NULL ) {
                              ISAME( 10 ) = LCE( BS, BB, LBB );
                           } else {
                             ISAME( 10 ) = LCERES( 'ge', ' ', M, N, BS, BB, LDB );
                           }
                           ISAME( 11 ) = LDBS == LDB;

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 50
                              SAME = SAME && ISAME( I );
                              IF( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                           } // 50
                           if ( !SAME ) {
                              FATAL = true;
                              GO TO 150;
                           }

                           if ( !NULL ) {
                              if ( SNAME( 10: 11 ) == 'mm' ) {

                                 // Check the result.

                                 if ( LEFT ) {
                                   cmmch(TRANSA, 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
                                 } else {
                                    cmmch('N', TRANSA, M, N, N, ALPHA, B, NMAX, A, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
                                 }
                              } else if ( SNAME( 10: 11 ) == 'sm' ) {

                                 // Compute approximation to original
                                 // matrix.

                                 for (J = 1; J <= N; J++) { // 70
                                    for (I = 1; I <= M; I++) { // 60
                                       C( I, J ) = BB( I + ( J - 1 )* LDB )                                        BB( I + ( J - 1 )*LDB ) = ALPHA* B( I, J );
                                    } // 60
                                 } // 70

                                 if ( LEFT ) {
                                    cmmch(TRANSA, 'N', M, N, M, ONE, A, NMAX, C, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
                                 } else {
                                    cmmch('N', TRANSA, M, N, N, ONE, C, NMAX, A, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
                                 }
                              }
                              ERRMAX = MAX( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 150;
                           }

                        } // 80

                     } // 90

                  } // 100

               } // 110

            } // 120

         } // 130

      } // 140

      // Report result.

      if ( ERRMAX < THRESH ) {
         if (IORDER == 0) WRITE( NOUT, FMT = 10000 )SNAME, NC;
         if (IORDER == 1) WRITE( NOUT, FMT = 10001 )SNAME, NC;
      } else {
         if (IORDER == 0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX;
         if (IORDER == 1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX;
      }
      GO TO 160;

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      if (TRACE) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);

      } // 160
      return;

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9996 FORMAT(' ******* ', A12,' FAILED ON CALL NUMBER:' );
 9995 FORMAT(1X, I6, ': ', A12,'(', 4( '''', A1, ''',' ), 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ')         ', '      .' );
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CCHK3.

      }

      SUBROUTINE CPRCN3(NOUT, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);
      int              NOUT, NC, IORDER, M, N, LDA, LDB;
      COMPLEX          ALPHA;
      String           SIDE, UPLO, TRANSA, DIAG;
      String           SNAME;
      String           CRC, CS, CU, CA, CD;

      if (SIDE == 'L') {
         CS =  '     CblasLeft';
      } else {
         CS =  '    CblasRight';
      }
      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (TRANSA == 'N') {
         CA =  '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CA =  '    CblasTrans';
      } else {
         CA =  'CblasConjTrans';
      }
      if (DIAG == 'N') {
         CD =  '  CblasNonUnit';
      } else {
         CD =  '     CblasUnit';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CS,CU;
      WRITE(NOUT, FMT = 9994)CA, CD, M, N, ALPHA, LDA, LDB;

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',');
 9994 FORMAT( 10X, 2( A14, ',') , 2( I3, ',' ), ' (', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ').' );
      }

      SUBROUTINE CCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER );

// Tests CHERK and CSYRK.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RONE, RZERO;
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS;
      REAL               ERR, ERRMAX, RALPHA, RALS, RBETA, RBETS;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, KS, LAA, LCC, LDA, LDAS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               CONJ, NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, TRANST, UPLO, UPLOS;
      String             ICHT, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CCHERK, CMAKE, CMMCH, CCSYRK
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NC'/, ICHU/'UL'/;
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ) == 'he';

      NARGS = 10;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN );
         // Set LDC to 1 more than minimum value if room.
         LDC = N;
         if (LDC < NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC > NMAX) GO TO 100;
         LCC = LDC*N;

         for (IK = 1; IK <= NIDIM; IK++) { // 90
            K = IDIM( IK );

            for (ICT = 1; ICT <= 2; ICT++) { // 80
               TRANS = ICHT( ICT: ICT );
               TRAN = TRANS == 'C';
               if (TRAN && !CONJ) TRANS = 'T';
               if ( TRAN ) {
                  MA = K;
                  NA = N;
               } else {
                  MA = N;
                  NA = K;
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA;
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 80;
               LAA = LDA*NA;

               // Generate the matrix A.

               cmake('ge', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );
                     if ( CONJ ) {
                        RALPHA = REAL( ALPHA );
                        ALPHA = CMPLX( RALPHA, RZERO );
                     }

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );
                        if ( CONJ ) {
                           RBETA = REAL( BETA );
                           BETA = CMPLX( RBETA, RZERO );
                        }
                        NULL = N <= 0;
                        if (CONJ) NULL = NULL || ( ( K <= 0 || RALPHA == RZERO ) && RBETA == RONE );

                        // Generate the matrix C.

                        cmake(SNAME( 8: 9 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1;

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO;
                        TRANSS = TRANS;
                        NS = N;
                        KS = K;
                        if ( CONJ ) {
                           RALS = RALPHA;
                        } else {
                           ALS = ALPHA;
                        }
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I );
                        } // 10
                        LDAS = LDA;
                        if ( CONJ ) {
                           RBETS = RBETA;
                        } else {
                           BETS = BETA;
                        }
                        for (I = 1; I <= LCC; I++) { // 20
                           CS( I ) = CC( I );
                        } // 20
                        LDCS = LDC;

                        // Call the subroutine.

                        if ( CONJ ) {
                           if (TRACE) CALL CPRCN6( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC);
                           if (REWI) REWIND NTRA;
                           ccherk(IORDER, UPLO, TRANS, N, K, RALPHA, AA, LDA, RBETA, CC, LDC );
                        } else {
                           if (TRACE) CALL CPRCN4( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
                           if (REWI) REWIND NTRA;
                           ccsyrk(IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 120;
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS == UPLO;
                        ISAME( 2 ) = TRANSS == TRANS;
                        ISAME( 3 ) = NS == N;
                        ISAME( 4 ) = KS == K;
                        if ( CONJ ) {
                           ISAME( 5 ) = RALS == RALPHA;
                        } else {
                           ISAME( 5 ) = ALS == ALPHA;
                        }
                        ISAME( 6 ) = LCE( AS, AA, LAA );
                        ISAME( 7 ) = LDAS == LDA;
                        if ( CONJ ) {
                           ISAME( 8 ) = RBETS == RBETA;
                        } else {
                           ISAME( 8 ) = BETS == BETA;
                        }
                        if ( NULL ) {
                           ISAME( 9 ) = LCE( CS, CC, LCC );
                        } else {
                           ISAME( 9 ) = LCERES( SNAME( 8: 9 ), UPLO, N, N, CS, CC, LDC );
                        }
                        ISAME( 10 ) = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 30
                           SAME = SAME && ISAME( I );
                           IF( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                        } // 30
                        if ( !SAME ) {
                           FATAL = true;
                           GO TO 120;
                        }

                        if ( !NULL ) {

                           // Check the result column by column.

                           if ( CONJ ) {
                              TRANST = 'C';
                           } else {
                              TRANST = 'T';
                           }
                           JC = 1;
                           for (J = 1; J <= N; J++) { // 40
                              if ( UPPER ) {
                                 JJ = 1;
                                 LJ = J;
                              } else {
                                 JJ = J;
                                 LJ = N - J + 1;
                              }
                              if ( TRAN ) {
                                 cmmch(TRANST, 'N', LJ, 1, K, ALPHA, A( 1, JJ ), NMAX, A( 1, J ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              } else {
                                 cmmch('N', TRANST, LJ, 1, K, ALPHA, A( JJ, 1 ), NMAX, A( J, 1 ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC;
                              } else {
                                 JC = JC + LDC + 1;
                              }
                              ERRMAX = MAX( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 110;
                           } // 40
                        }

                     } // 50

                  } // 60

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX < THRESH ) {
         if (IORDER == 0) WRITE( NOUT, FMT = 10000 )SNAME, NC;
         if (IORDER == 1) WRITE( NOUT, FMT = 10001 )SNAME, NC;
      } else {
         if (IORDER == 0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX;
         if (IORDER == 1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 110
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( CONJ ) {
      cprcn6(NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, RALPHA, LDA, rBETA, LDC);
      } else {
      cprcn4(NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
      }

      } // 130
      return;

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );
 9994 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ), F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ')               ', '          .' );
 9993 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ), '(', F4.1, ',', F4.1, ') , A,', I3, ',(', F4.1, ',', F4.1, '), C,', I3, ')          .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CCHK4.

      }

      SUBROUTINE CPRCN4(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, BETA, LDC);
      int              NOUT, NC, IORDER, N, K, LDA, LDC;
      COMPLEX          ALPHA, BETA;
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (TRANSA == 'N') {
         CA =  '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CA =  '    CblasTrans';
      } else {
         CA =  'CblasConjTrans';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA;
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') );
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1 ,'), A,', I3, ', (', F4.1,',', F4.1, '), C,', I3, ').' );
      }


      SUBROUTINE CPRCN6(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, BETA, LDC);
      int              NOUT, NC, IORDER, N, K, LDA, LDC;
      REAL             ALPHA, BETA;
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (TRANSA == 'N') {
         CA =  '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CA =  '    CblasTrans';
      } else {
         CA =  'CblasConjTrans';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA;
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') );
 9994 FORMAT( 10X, 2( I3, ',' ), F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ').' );
      }

      SUBROUTINE CCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, IORDER );

// Tests CHER2K and CSYR2K.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RONE, RZERO;
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX );
      REAL               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS;
      REAL               ERR, ERRMAX, RBETA, RBETS;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               CONJ, NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, TRANST, UPLO, UPLOS;
      String             ICHT, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CCHER2K, CMAKE, CMMCH, CCSYR2K
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, MAX, REAL
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NC'/, ICHU/'UL'/;
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ) == 'he';

      NARGS = 12;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 130
         N = IDIM( IN );
         // Set LDC to 1 more than minimum value if room.
         LDC = N;
         if (LDC < NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC > NMAX) GO TO 130;
         LCC = LDC*N;

         for (IK = 1; IK <= NIDIM; IK++) { // 120
            K = IDIM( IK );

            for (ICT = 1; ICT <= 2; ICT++) { // 110
               TRANS = ICHT( ICT: ICT );
               TRAN = TRANS == 'C';
               if (TRAN && !CONJ) TRANS = 'T';
               if ( TRAN ) {
                  MA = K;
                  NA = N;
               } else {
                  MA = N;
                  NA = K;
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA;
               if (LDA < NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA > NMAX) GO TO 110;
               LAA = LDA*NA;

               // Generate the matrix A.

               if ( TRAN ) {
                  cmake('ge', ' ', ' ', MA, NA, AB, 2*NMAX, AA, LDA, RESET, ZERO );
               } else {
                 cmake('ge', ' ', ' ', MA, NA, AB, NMAX, AA, LDA, RESET, ZERO );
               }

               // Generate the matrix B.

               LDB = LDA;
               LBB = LAA;
               if ( TRAN ) {
                  cmake('ge', ' ', ' ', MA, NA, AB( K + 1 ), 2*NMAX, BB, LDB, RESET, ZERO );
               } else {
                  cmake('ge', ' ', ' ', MA, NA, AB( K*NMAX + 1 ), NMAX, BB, LDB, RESET, ZERO );
               }

               for (ICU = 1; ICU <= 2; ICU++) { // 100
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 90
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 80
                        BETA = BET( IB );
                        if ( CONJ ) {
                           RBETA = REAL( BETA );
                           BETA = CMPLX( RBETA, RZERO );
                        }
                        NULL = N <= 0;
                        if (CONJ) NULL = NULL || ( ( K <= 0 || ALPHA == ZERO ) && RBETA == RONE );

                        // Generate the matrix C.

                        cmake(SNAME( 8: 9 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1;

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO;
                        TRANSS = TRANS;
                        NS = N;
                        KS = K;
                        ALS = ALPHA;
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I );
                        } // 10
                        LDAS = LDA;
                        for (I = 1; I <= LBB; I++) { // 20
                           BS( I ) = BB( I );
                        } // 20
                        LDBS = LDB;
                        if ( CONJ ) {
                           RBETS = RBETA;
                        } else {
                           BETS = BETA;
                        }
                        for (I = 1; I <= LCC; I++) { // 30
                           CS( I ) = CC( I );
                        } // 30
                        LDCS = LDC;

                        // Call the subroutine.

                        if ( CONJ ) {
                           if (TRACE) CALL CPRCN7( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC);
                           if (REWI) REWIND NTRA;
                           ccher2k(IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, RBETA, CC, LDC );
                        } else {
                           if (TRACE) CALL CPRCN5( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
                           if (REWI) REWIND NTRA;
                           ccsyr2k(IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 150;
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS == UPLO;
                        ISAME( 2 ) = TRANSS == TRANS;
                        ISAME( 3 ) = NS == N;
                        ISAME( 4 ) = KS == K;
                        ISAME( 5 ) = ALS == ALPHA;
                        ISAME( 6 ) = LCE( AS, AA, LAA );
                        ISAME( 7 ) = LDAS == LDA;
                        ISAME( 8 ) = LCE( BS, BB, LBB );
                        ISAME( 9 ) = LDBS == LDB;
                        if ( CONJ ) {
                           ISAME( 10 ) = RBETS == RBETA;
                        } else {
                           ISAME( 10 ) = BETS == BETA;
                        }
                        if ( NULL ) {
                           ISAME( 11 ) = LCE( CS, CC, LCC );
                        } else {
                           ISAME( 11 ) = LCERES( 'he', UPLO, N, N, CS, CC, LDC );
                        }
                        ISAME( 12 ) = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           IF( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                        } // 40
                        if ( !SAME ) {
                           FATAL = true;
                           GO TO 150;
                        }

                        if ( !NULL ) {

                           // Check the result column by column.

                           if ( CONJ ) {
                              TRANST = 'C';
                           } else {
                              TRANST = 'T';
                           }
                           JJAB = 1;
                           JC = 1;
                           for (J = 1; J <= N; J++) { // 70
                              if ( UPPER ) {
                                 JJ = 1;
                                 LJ = J;
                              } else {
                                 JJ = J;
                                 LJ = N - J + 1;
                              }
                              if ( TRAN ) {
                                 for (I = 1; I <= K; I++) { // 50
                                    W( I ) = ALPHA*AB( ( J - 1 )*2* NMAX + K + I );
                                    if ( CONJ ) {
                                       W( K + I ) = CONJG( ALPHA )* AB( ( J - 1 )*2* NMAX + I );
                                    } else {
                                       W( K + I ) = ALPHA* AB( ( J - 1 )*2* NMAX + I );
                                    }
                                 } // 50
                                 cmmch(TRANST, 'N', LJ, 1, 2*K, ONE, AB( JJAB ), 2*NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              } else {
                                 for (I = 1; I <= K; I++) { // 60
                                    if ( CONJ ) {
                                       W( I ) = ALPHA*CONJG( AB( ( K + I - 1 )*NMAX + J ) )                                        W( K + I ) = CONJG( ALPHA* AB( ( I - 1 )*NMAX + J ) );
                                    } else {
                                       W( I ) = ALPHA*AB( ( K + I - 1 )* NMAX + J )                                        W( K + I ) = ALPHA* AB( ( I - 1 )*NMAX + J );
                                    }
                                 } // 60
                                 cmmch('N', 'N', LJ, 1, 2*K, ONE, AB( JJ ), NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC;
                              } else {
                                 JC = JC + LDC + 1;
                                 if (TRAN) JJAB = JJAB + 2*NMAX;
                              }
                              ERRMAX = MAX( ERRMAX, ERR );
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 140;
                           } // 70
                        }

                     } // 80

                  } // 90

               } // 100

            } // 110

         } // 120

      } // 130

      // Report result.

      if ( ERRMAX < THRESH ) {
         if (IORDER == 0) WRITE( NOUT, FMT = 10000 )SNAME, NC;
         if (IORDER == 1) WRITE( NOUT, FMT = 10001 )SNAME, NC;
      } else {
         if (IORDER == 0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX;
         if (IORDER == 1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX;
      }
      GO TO 160;

      } // 140
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( CONJ ) {
         cprcn7(NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC);
      } else {
         cprcn5(NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
      }

      } // 160
      return;

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ', 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ', 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS', ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' );
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );
 9994 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',', F4.1, ', C,', I3, ')           .' );
 9993 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1, ',', F4.1, '), C,', I3, ')    .' );
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' );

      // End of CCHK5.

      }

      SUBROUTINE CPRCN5(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, LDB, BETA, LDC);
      int              NOUT, NC, IORDER, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA;
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (TRANSA == 'N') {
         CA =  '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CA =  '    CblasTrans';
      } else {
         CA =  'CblasConjTrans';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA;
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, LDB, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') );
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1, '), A,', I3, ', B', I3, ', (', F4.1, ',', F4.1, '), C,', I3, ').' );
      }


      SUBROUTINE CPRCN7(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, LDB, BETA, LDC);
      int              NOUT, NC, IORDER, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA;
      REAL             BETA;
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO == 'U') {
         CU =  '    CblasUpper';
      } else {
         CU =  '    CblasLower';
      }
      if (TRANSA == 'N') {
         CA =  '  CblasNoTrans';
      } else if (TRANSA == 'T') {
         CA =  '    CblasTrans';
      } else {
         CA =  'CblasConjTrans';
      }
      if (IORDER == 1) {
         CRC = ' CblasRowMajor';
      } else {
         CRC = ' CblasColMajor';
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA;
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, LDB, BETA, LDC;

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') );
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1, '), A,', I3, ', B', I3, ',', F4.1, ', C,', I3, ').' );
      }

      SUBROUTINE CMAKE(TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, RESET, TRANSL );

// Generates values for an M by N matrix A.
// Stores the values in the array AA in the data structure required
// by the routine, with unwanted elements set to rogue value.

// TYPE is 'ge', 'he', 'sy' or 'tr'.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      COMPLEX            ROGUE;
      const              ROGUE = ( -1.0e10, 1.0e10 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      REAL               RROGUE;
      const              RROGUE = -1.0e10 ;
      // .. Scalar Arguments ..
      COMPLEX            TRANSL;
      int                LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J, JJ;
      bool               GEN, HER, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      COMPLEX            CBEG;
      // EXTERNAL CBEG
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, REAL
      // .. Executable Statements ..
      GEN = TYPE == 'ge';
      HER = TYPE == 'he';
      SYM = TYPE == 'sy';
      TRI = TYPE == 'tr';
      UPPER = ( HER || SYM || TRI ) && UPLO == 'U';
      LOWER = ( HER || SYM || TRI ) && UPLO == 'L';
      UNIT = TRI && DIAG == 'U';

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN || ( UPPER && I <= J ) || ( LOWER && I >= J ) ) {
               A( I, J ) = CBEG( RESET ) + TRANSL;
               if ( I != J ) {
                  // Set some elements to zero
                  if (N > 3 && J == N/2) A( I, J ) = ZERO;
                  if ( HER ) {
                     A( J, I ) = CONJG( A( I, J ) );
                  } else if ( SYM ) {
                     A( J, I ) = A( I, J );
                  } else if ( TRI ) {
                     A( J, I ) = ZERO;
                  }
               }
            }
         } // 10
         if (HER) A( J, J ) = CMPLX( REAL( A( J, J ) ), RZERO )          IF( TRI ) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE;
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
      } else if ( TYPE == 'he' || TYPE == 'sy' || TYPE == 'tr' ) {
         for (J = 1; J <= N; J++) { // 90
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
            for (I = 1; I <= IBEG - 1; I++) { // 60
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 60
            for (I = IBEG; I <= IEND; I++) { // 70
               AA( I + ( J - 1 )*LDA ) = A( I, J );
            } // 70
            for (I = IEND + 1; I <= LDA; I++) { // 80
               AA( I + ( J - 1 )*LDA ) = ROGUE;
            } // 80
            if ( HER ) {
               JJ = J + ( J - 1 )*LDA;
               AA( JJ ) = CMPLX( REAL( AA( JJ ) ), RROGUE );
            }
         } // 90
      }
      return;

      // End of CMAKE.

      }
      SUBROUTINE CMMCH(TRANSA, TRANSB, M, N, KK, ALPHA, A, LDA, B, LDB, BETA, C, LDC, CT, G, CC, LDCC, EPS, ERR, FATAL, NOUT, MV );

// Checks the results of the computational tests.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO, RONE;
      const              RZERO = 0.0, RONE = 1.0 ;
      // .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA;
      REAL               EPS, ERR;
      int                KK, LDA, LDB, LDC, LDCC, M, N, NOUT;
      bool               FATAL, MV;
      String             TRANSA, TRANSB;
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), CC( LDCC, * ), CT( * );
      REAL               G( * );
      // .. Local Scalars ..
      COMPLEX            CL;
      REAL               ERRI;
      int                I, J, K;
      bool               CTRANA, CTRANB, TRANA, TRANB;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
      // .. Statement Functions ..
      REAL               ABS1;
      // .. Statement Function definitions ..
      ABS1( CL ) = ABS( REAL( CL ) ) + ABS( AIMAG( CL ) );
      // .. Executable Statements ..
      TRANA = TRANSA == 'T' || TRANSA == 'C';
      TRANB = TRANSB == 'T' || TRANSB == 'C';
      CTRANA = TRANSA == 'C';
      CTRANB = TRANSB == 'C';

      // Compute expected result, one column at a time, in CT using data
      // in A, B and C.
      // Compute gauges in G.

      for (J = 1; J <= N; J++) { // 220

         for (I = 1; I <= M; I++) { // 10
            CT( I ) = ZERO;
            G( I ) = RZERO;
         } // 10
         if ( !TRANA && !TRANB ) {
            for (K = 1; K <= KK; K++) { // 30
               for (I = 1; I <= M; I++) { // 20
                  CT( I ) = CT( I ) + A( I, K )*B( K, J );
                  G( I ) = G( I ) + ABS1( A( I, K ) )*ABS1( B( K, J ) );
               } // 20
            } // 30
         } else if ( TRANA && !TRANB ) {
            if ( CTRANA ) {
               for (K = 1; K <= KK; K++) { // 50
                  for (I = 1; I <= M; I++) { // 40
                     CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( K, J );
                     G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) );
                  } // 40
               } // 50
            } else {
               for (K = 1; K <= KK; K++) { // 70
                  for (I = 1; I <= M; I++) { // 60
                     CT( I ) = CT( I ) + A( K, I )*B( K, J );
                     G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) );
                  } // 60
               } // 70
            }
         } else if ( !TRANA && TRANB ) {
            if ( CTRANB ) {
               for (K = 1; K <= KK; K++) { // 90
                  for (I = 1; I <= M; I++) { // 80
                     CT( I ) = CT( I ) + A( I, K )*CONJG( B( J, K ) );
                     G( I ) = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) );
                  } // 80
               } // 90
            } else {
               for (K = 1; K <= KK; K++) { // 110
                  for (I = 1; I <= M; I++) { // 100
                     CT( I ) = CT( I ) + A( I, K )*B( J, K );
                     G( I ) = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) );
                  } // 100
               } // 110
            }
         } else if ( TRANA && TRANB ) {
            if ( CTRANA ) {
               if ( CTRANB ) {
                  for (K = 1; K <= KK; K++) { // 130
                     for (I = 1; I <= M; I++) { // 120
                        CT( I ) = CT( I ) + CONJG( A( K, I ) )* CONJG( B( J, K ) )                         G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 120
                  } // 130
               } else {
                  for (K = 1; K <= KK; K++) { // 150
                     for (I = 1; I <= M; I++) { // 140
                       CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( J, K );
                       G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 140
                  } // 150
               }
            } else {
               if ( CTRANB ) {
                  for (K = 1; K <= KK; K++) { // 170
                     for (I = 1; I <= M; I++) { // 160
                       CT( I ) = CT( I ) + A( K, I )*CONJG( B( J, K ) );
                       G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 160
                  } // 170
               } else {
                  for (K = 1; K <= KK; K++) { // 190
                     for (I = 1; I <= M; I++) { // 180
                        CT( I ) = CT( I ) + A( K, I )*B( J, K );
                        G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 180
                  } // 190
               }
            }
         }
         for (I = 1; I <= M; I++) { // 200
            CT( I ) = ALPHA*CT( I ) + BETA*C( I, J );
            G( I ) = ABS1( ALPHA )*G( I ) + ABS1( BETA )*ABS1( C( I, J ) );
         } // 200

         // Compute the error ratio for this result.

         ERR = ZERO;
         for (I = 1; I <= M; I++) { // 210
            ERRI = ABS1( CT( I ) - CC( I, J ) )/EPS;
            IF( G( I ) != RZERO ) ERRI = ERRI/G( I );
            ERR = MAX( ERR, ERRI );
            IF( ERR*SQRT( EPS ) >= RONE ) GO TO 230;
         } // 210

      } // 220

      // If the loop completes, all results are at least half accurate.
      GO TO 250;

      // Report fatal error.

  230 FATAL = true;
      WRITE( NOUT, FMT = 9999 );
      for (I = 1; I <= M; I++) { // 240
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, CT( I ), CC( I, J );
         } else {
            WRITE( NOUT, FMT = 9998 )I, CC( I, J ), CT( I );
         }
      } // 240
      if (N > 1) WRITE( NOUT, FMT = 9997 )J;

      } // 250
      return;

 9999 FORMAT(' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'                       EXPECTED RE', 'SULT                    COMPUTED RESULT' );
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) );
 9997 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 );

      // End of CMMCH.

      }
      bool    FUNCTION LCE( RI, RJ, LR );

// Tests if two arrays are identical.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      COMPLEX            RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         IF( RI( I ) != RJ( I ) ) GO TO 20;
      } // 10
      LCE = true;
      GO TO 30;
      } // 20
      LCE = false;
   30 RETURN;

      // End of LCE.

      }
      bool    FUNCTION LCERES( TYPE, UPLO, M, N, AA, AS, LDA );

// Tests if selected elements in two arrays are equal.

// TYPE is 'ge' or 'he' or 'sy'.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX            AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO == 'U';
      if ( TYPE == 'ge' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 10
         } // 20
      } else if ( TYPE == 'he' || TYPE == 'sy' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1;
               IEND = J;
            } else {
               IBEG = J;
               IEND = N;
            }
            for (I = 1; I <= IBEG - 1; I++) { // 30
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 30
            for (I = IEND + 1; I <= LDA; I++) { // 40
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 40
         } // 50
      }

      } // 60
      LCERES = true;
      GO TO 80;
      } // 70
      LCERES = false;
   80 RETURN;

      // End of LCERES.

      }
      COMPLEX FUNCTION CBEG( RESET );

// Generates complex numbers as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ;
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX
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
      CBEG = CMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 );
      return;

      // End of CBEG.

      }
      REAL FUNCTION SDIFF( X, Y );

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      REAL               X, Y;
      // .. Executable Statements ..
      SDIFF = X - Y;
      return;

      // End of SDIFF.

      }
