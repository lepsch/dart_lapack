      void main() {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 9 ;
      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      int                NMAX;
      const              NMAX = 65 ;
      int                NIDMAX, NALMAX, NBEMAX;
      const              NIDMAX = 9, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      double               EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NOUT, NTRA;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR;
      String             TRANSA, TRANSB;
      String             SNAMET;
      String             SNAPS, SUMMRY;
      // .. Local Arrays ..
      Complex            AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBEMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX );
      double               G( NMAX );
      int                IDIM( NIDMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      //- REAL               SDIFF;
      //- bool               LCE;
      // EXTERNAL SDIFF, LCE
      // .. External Subroutines ..
      // EXTERNAL CCHK1, CCHK2, CCHK3, CCHK4, CCHK5, CCHKE, CMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // COMMON /SRNAMC(srnamc.SRNAMT
      // .. Data statements ..
      const SNAMES = ['CGEMM ', 'CHEMM ', 'CSYMM ', 'CTRMM ', 'CTRSM ', 'CHERK ', 'CSYRK ', 'CHER2K', 'CSYR2K'];


      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SUMMRY;
      READ( NIN, FMT = * )NOUT;
      OPEN( NOUT, FILE = SUMMRY );
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

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      for (I = 1; I <= NSUBS; I++) { // 20
         LTEST[I] = false;
      } // 20
   30 READ( NIN, FMT = 9988, END = 60 )SNAMET, LTESTT;
      for (I = 1; I <= NSUBS; I++) { // 40
         if( SNAMET == SNAMES( I ) ) GO TO 50;
      } // 40
      WRITE( NOUT, FMT = 9990 )SNAMET;
      STOP;
   50 LTEST( I ) = LTESTT;
      GO TO 30;

      } // 60
      CLOSE ( NIN );

      // Compute EPS (the machine precision).

      EPS = EPSILON(RZERO);
      WRITE( NOUT, FMT = 9998 )EPS;

      // Check the reliability of CMMCH using exact data.

      N = min( 32, NMAX );
      for (J = 1; J <= N; J++) { // 100
         for (I = 1; I <= N; I++) { // 90
            AB[I][J] = max( I - J + 1, 0 );
         } // 90
         AB[J][NMAX + 1] = J;
         AB[1][NMAX + J] = J;
         C[J][1] = ZERO;
      } // 100
      for (J = 1; J <= N; J++) { // 110
         CC[J] = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
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
         AB[J][NMAX + 1] = N - J + 1;
         AB[1][NMAX + J] = N - J + 1;
      } // 120
      for (J = 1; J <= N; J++) { // 130
         CC[N - J + 1] = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3;
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
           srnamc.SRNAMT = SNAMES( ISNUM );
            // Test error exits.
            if ( TSTERR ) {
               cchke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * );
            }
            // Test computations.
            INFOT = 0;
            OK = true;
            FATAL = false;
            GO TO ( 140, 150, 150, 160, 160, 170, 170, 180, 180 )ISNUM;
            // Test CGEMM, 01.
  140       CALL CCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test CHEMM, 02, CSYMM, 03.
  150       CALL CCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test CTRMM, 04, CTRSM, 05.
  160       CALL CCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C );
            GO TO 190;
            // Test CHERK, 06, CSYRK, 07.
  170       CALL CCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test CHER2K, 08, CSYR2K, 09.
  180       CALL CCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W );
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

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LESS THAN${.f8_2}');
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 );
 9997 FORMAT( ' NUMBER OF VALUES OF ${} IS LESS THAN 1 OR GREATER THAN ${.i2}');
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${.i2}');
 9995 FORMAT( ' TESTS OF THE COMPLEX          LEVEL 3 BLAS', //' THE FOLLOWING PARAMETER VALUES WILL BE USED:' );
 9994 FORMAT( '   FOR N              ', 9I6 );
 9993 FORMAT( '   FOR ALPHA          ', 7( '(${.f4_1},${.f4_1})  ', : ) );
 9992 FORMAT( '   FOR BETA           ', 7( '(${.f4_1},${.f4_1})  ', : ) );
 9991 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******' );
 9990 FORMAT( ' SUBPROGRAM NAME ${.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******' );
 9989 FORMAT( ' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n CMMCH WAS CALLED WITH TRANSA = ${.a1} AND TRANSB = ${.a1}\n AND RETURNED SAME = ${.l1} AND ERR = ${.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******' );
 9988 FORMAT( A6, L2 );
 9987 FORMAT(' ${.a6} WAS NOT TESTED' );
 9986 FORMAT('\n END OF TESTS' );
 9985 FORMAT('\n ******* FATAL ERROR - TESTS ABANDONED *******' );
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' );
      }
      void cchk1(final int SNAME, EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NBET, final int BET, final int NMAX, final int A, final int AA, final int AS, final int B, final int BB, final int BS, final int C, final int CC, final int CS, final int CT, final int G,) {

// Tests CGEMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      double               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      Complex            ALPHA, ALS, BETA, BLS;
      double               ERR, ERRMAX;
      int                I, IA, IB, ICA, ICB, IK, IM, IN, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MA, MB, MS, N, NA, NARGS, NB, NC, NS;
      bool               NULL, RESET, SAME, TRANA, TRANB;
      String             TRANAS, TRANBS, TRANSA, TRANSB;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CMAKE, CMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICH = 'NTC';


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

                  cmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

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

                     cmake('GE', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO );

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA );

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB );

                           // Generate the matrix C.

                           cmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                              AS[I] = AA( I );
                           } // 10
                           LDAS = LDA;
                           for (I = 1; I <= LBB; I++) { // 20
                              BS[I] = BB( I );
                           } // 20
                           LDBS = LDB;
                           BLS = BETA;
                           for (I = 1; I <= LCC; I++) { // 30
                              CS[I] = CC( I );
                           } // 30
                           LDCS = LDC;

                           // Call the subroutine.

                           if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC;
                           if (REWI) REWIND NTRA;
                           cgemm(TRANSA, TRANSB, M, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                           // Check if error-exit was taken incorrectly.

                           if ( !OK ) {
                              WRITE( NOUT, FMT = 9994 );
                              FATAL = true;
                              GO TO 120;
                           }

                           // See what data changed inside subroutines.

                           ISAME[1] = TRANSA == TRANAS;
                           ISAME[2] = TRANSB == TRANBS;
                           ISAME[3] = MS == M;
                           ISAME[4] = NS == N;
                           ISAME[5] = KS == K;
                           ISAME[6] = ALS == ALPHA;
                           ISAME[7] = LCE( AS, AA, LAA );
                           ISAME[8] = LDAS == LDA;
                           ISAME[9] = LCE( BS, BB, LBB );
                           ISAME[10] = LDBS == LDB;
                           ISAME[11] = BLS == BETA;
                           if ( NULL ) {
                              ISAME[12] = LCE( CS, CC, LCC );
                           } else {
                              ISAME[12] = LCERES( 'GE', ' ', M, N, CS, CC, LDC );
                           }
                           ISAME[13] = LDCS == LDC;

                           // If data was incorrectly changed, report
                           // and return.

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

                              cmmch(TRANSA, TRANSB, M, N, K, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
                              ERRMAX = max( ERRMAX, ERR );
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
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC;

      } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(\'${.a1}\',\'${.a1}\',${i3(3, ',')}', '(${.f4_1},${.f4_1}), A,${.i3}, B,${.i3},(${.f4_1},${.f4_1}), C,${.i3}).' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void cchk2(final int SNAME, EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NBET, final int BET, final int NMAX, final int A, final int AA, final int AS, final int B, final int BB, final int BS, final int C, final int CC, final int CS, final int CT, final int G,) {

// Tests CHEMM and CSYMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      double               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      Complex            ALPHA, ALS, BETA, BLS;
      double               ERR, ERRMAX;
      int                I, IA, IB, ICS, ICU, IM, IN, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MS, N, NA, NARGS, NC, NS;
      bool               CONJ, LEFT, NULL, RESET, SAME;
      String             SIDE, SIDES, UPLO, UPLOS;
      String             ICHS, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CHEMM, CMAKE, CMMCH, CSYMM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHS = 'LR', ICHU = 'UL';

      CONJ = SNAME( 2: 3 ) == 'HE';

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

            cmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

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

                  cmake(SNAME( 2: 3 ), UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );

                        // Generate the matrix C.

                        cmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1;

                        // Save every datum before calling the
                        // subroutine.

                        SIDES = SIDE;
                        UPLOS = UPLO;
                        MS = M;
                        NS = N;
                        ALS = ALPHA;
                        for (I = 1; I <= LAA; I++) { // 10
                           AS[I] = AA( I );
                        } // 10
                        LDAS = LDA;
                        for (I = 1; I <= LBB; I++) { // 20
                           BS[I] = BB( I );
                        } // 20
                        LDBS = LDB;
                        BLS = BETA;
                        for (I = 1; I <= LCC; I++) { // 30
                           CS[I] = CC( I );
                        } // 30
                        LDCS = LDC;

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        if ( CONJ ) {
                           chemm(SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        } else {
                           csymm(SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9994 );
                           FATAL = true;
                           GO TO 110;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = SIDES == SIDE;
                        ISAME[2] = UPLOS == UPLO;
                        ISAME[3] = MS == M;
                        ISAME[4] = NS == N;
                        ISAME[5] = ALS == ALPHA;
                        ISAME[6] = LCE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        ISAME[8] = LCE( BS, BB, LBB );
                        ISAME[9] = LDBS == LDB;
                        ISAME[10] = BLS == BETA;
                        if ( NULL ) {
                           ISAME[11] = LCE( CS, CC, LCC );
                        } else {
                           ISAME[11] = LCERES( 'GE', ' ', M, N, CS, CC, LDC );
                        }
                        ISAME[12] = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
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
                           ERRMAX = max( ERRMAX, ERR );
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
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 120;

      } // 110
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC;

      } // 120
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), '(${.f4_1},${.f4_1}), A,${.i3}, B,${.i3},(${.f4_1},${.f4_1}), C,${.i3})    .' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void cchk3(final int SNAME, final int EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NMAX, final int A, final int AA, final int AS, final int B, final int BB, final int BS, final int CT, final int G, final int C,) {

// Tests CTRMM and CTRSM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CT( NMAX );
      double               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      Complex            ALPHA, ALS;
      double               ERR, ERRMAX;
      int                I, IA, ICD, ICS, ICT, ICU, IM, IN, J, LAA, LBB, LDA, LDAS, LDB, LDBS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String             DIAG, DIAGS, SIDE, SIDES, TRANAS, TRANSA, UPLO, UPLOS;
      String             ICHD, ICHS, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CMAKE, CMMCH, CTRMM, CTRSM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN', ICHS = 'LR';


      NARGS = 11;
      NC = 0;
      RESET = true;
      ERRMAX = RZERO;
      // Set up zero matrix for CMMCH.
      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            C[I][J] = ZERO;
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

                           cmake('TR', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                           // Generate the matrix B.

                           cmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

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
                              AS[I] = AA( I );
                           } // 30
                           LDAS = LDA;
                           for (I = 1; I <= LBB; I++) { // 40
                              BS[I] = BB( I );
                           } // 40
                           LDBS = LDB;

                           // Call the subroutine.

                           if ( SNAME( 4: 5 ) == 'MM' ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;
                              if (REWI) REWIND NTRA;
                              ctrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           } else if ( SNAME( 4: 5 ) == 'SM' ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;
                              if (REWI) REWIND NTRA;
                              ctrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( !OK ) {
                              WRITE( NOUT, FMT = 9994 );
                              FATAL = true;
                              GO TO 150;
                           }

                           // See what data changed inside subroutines.

                           ISAME[1] = SIDES == SIDE;
                           ISAME[2] = UPLOS == UPLO;
                           ISAME[3] = TRANAS == TRANSA;
                           ISAME[4] = DIAGS == DIAG;
                           ISAME[5] = MS == M;
                           ISAME[6] = NS == N;
                           ISAME[7] = ALS == ALPHA;
                           ISAME[8] = LCE( AS, AA, LAA );
                           ISAME[9] = LDAS == LDA;
                           if ( NULL ) {
                              ISAME[10] = LCE( BS, BB, LBB );
                           } else {
                              ISAME[10] = LCERES( 'GE', ' ', M, N, BS, BB, LDB );
                           }
                           ISAME[11] = LDBS == LDB;

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 50
                              SAME = SAME && ISAME( I );
                              if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
                           } // 50
                           if ( !SAME ) {
                              FATAL = true;
                              GO TO 150;
                           }

                           if ( !NULL ) {
                              if ( SNAME( 4: 5 ) == 'MM' ) {

                                 // Check the result.

                                 if ( LEFT ) {
                                    cmmch(TRANSA, 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
                                 } else {
                                    cmmch('N', TRANSA, M, N, N, ALPHA, B, NMAX, A, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
                                 }
                              } else if ( SNAME( 4: 5 ) == 'SM' ) {

                                 // Compute approximation to original
                                 // matrix.

                                 for (J = 1; J <= N; J++) { // 70
                                    for (I = 1; I <= M; I++) { // 60
                                       C[I][J] = BB( I + ( J - 1 )* LDB )                                        BB( I + ( J - 1 )*LDB ) = ALPHA* B( I, J );
                                    } // 60
                                 } // 70

                                 if ( LEFT ) {
                                    cmmch(TRANSA, 'N', M, N, M, ONE, A, NMAX, C, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
                                 } else {
                                    cmmch('N', TRANSA, M, N, N, ONE, C, NMAX, A, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
                                 }
                              }
                              ERRMAX = max( ERRMAX, ERR );
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
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 160;

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;

      } // 160
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(', 4( '\'${.a1}\',' ), 2( I3, ',' ), '(${.f4_1},${.f4_1}), A,${.i3}, B,${.i3})               .' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void cchk4(final int SNAME, EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NBET, final int BET, final int NMAX, final int A, final int AA, final int AS, final int B, final int BB, final int BS, final int C, final int CC, final int CS, final int CT, final int G,) {

// Tests CHERK and CSYRK.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double               RONE, RZERO;
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX );
      double               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      Complex            ALPHA, ALS, BETA, BETS;
      double               ERR, ERRMAX, RALPHA, RALS, RBETA, RBETS;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, KS, LAA, LCC, LDA, LDAS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               CONJ, NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, TRANST, UPLO, UPLOS;
      String             ICHT, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CHERK, CMAKE, CMMCH, CSYRK
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHT = 'NC', ICHU = 'UL';

      CONJ = SNAME( 2: 3 ) == 'HE';

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

               cmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );
                     if ( CONJ ) {
                        RALPHA = double( ALPHA );
                        ALPHA = CMPLX( RALPHA, RZERO );
                     }

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );
                        if ( CONJ ) {
                           RBETA = double( BETA );
                           BETA = CMPLX( RBETA, RZERO );
                        }
                        NULL = N <= 0;
                        if (CONJ) NULL = NULL || ( ( K <= 0 || RALPHA == RZERO ) && RBETA == RONE );

                        // Generate the matrix C.

                        cmake(SNAME( 2: 3 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                           AS[I] = AA( I );
                        } // 10
                        LDAS = LDA;
                        if ( CONJ ) {
                           RBETS = RBETA;
                        } else {
                           BETS = BETA;
                        }
                        for (I = 1; I <= LCC; I++) { // 20
                           CS[I] = CC( I );
                        } // 20
                        LDCS = LDC;

                        // Call the subroutine.

                        if ( CONJ ) {
                           if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC;
                           if (REWI) REWIND NTRA;
                           cherk(UPLO, TRANS, N, K, RALPHA, AA, LDA, RBETA, CC, LDC );
                        } else {
                           if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC;
                           if (REWI) REWIND NTRA;
                           csyrk(UPLO, TRANS, N, K, ALPHA, AA, LDA, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 120;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = UPLOS == UPLO;
                        ISAME[2] = TRANSS == TRANS;
                        ISAME[3] = NS == N;
                        ISAME[4] = KS == K;
                        if ( CONJ ) {
                           ISAME[5] = RALS == RALPHA;
                        } else {
                           ISAME[5] = ALS == ALPHA;
                        }
                        ISAME[6] = LCE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        if ( CONJ ) {
                           ISAME[8] = RBETS == RBETA;
                        } else {
                           ISAME[8] = BETS == BETA;
                        }
                        if ( NULL ) {
                           ISAME[9] = LCE( CS, CC, LCC );
                        } else {
                           ISAME[9] = LCERES( SNAME( 2: 3 ), UPLO, N, N, CS, CC, LDC );
                        }
                        ISAME[10] = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

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
                              ERRMAX = max( ERRMAX, ERR );
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
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 130;

      } // 110
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( CONJ ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC;
      } else {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC;
      }

      } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), F4.1, ', A,${.i3},${.f4_1}, C,${.i3})                         .' );
 9993 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), '(${.f4_1},${.f4_1}) , A,${.i3},(${.f4_1},${.f4_1}), C,${.i3})          .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void cchk5(final int SNAME, EPS, final int THRESH, final int NOUT, final int NTRA, final int TRACE, final int REWI, final int FATAL, final int NIDIM, final int IDIM, final int NALF, final int ALF, final int NBET, final int BET, final int NMAX, final int AB, final int AA, final int AS, final int BB, final int BS, final int C, final int CC, final int CS, final int CT, final int G, final int W,) {

// Tests CHER2K and CSYR2K.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double               RONE, RZERO;
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      Complex            AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX );
      double               G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      Complex            ALPHA, ALS, BETA, BETS;
      double               ERR, ERRMAX, RBETA, RBETS;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               CONJ, NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, TRANST, UPLO, UPLOS;
      String             ICHT, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LCE, LCERES;
      // EXTERNAL LCE, LCERES
      // .. External Subroutines ..
      // EXTERNAL CHER2K, CMAKE, CMMCH, CSYR2K
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, MAX, REAL
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHT = 'NC', ICHU = 'UL';

      CONJ = SNAME( 2: 3 ) == 'HE';

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
                  cmake('GE', ' ', ' ', MA, NA, AB, 2*NMAX, AA, LDA, RESET, ZERO );
               } else {
                  cmake('GE', ' ', ' ', MA, NA, AB, NMAX, AA, LDA, RESET, ZERO );
               }

               // Generate the matrix B.

               LDB = LDA;
               LBB = LAA;
               if ( TRAN ) {
                  cmake('GE', ' ', ' ', MA, NA, AB( K + 1 ), 2*NMAX, BB, LDB, RESET, ZERO );
               } else {
                  cmake('GE', ' ', ' ', MA, NA, AB( K*NMAX + 1 ), NMAX, BB, LDB, RESET, ZERO );
               }

               for (ICU = 1; ICU <= 2; ICU++) { // 100
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 90
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 80
                        BETA = BET( IB );
                        if ( CONJ ) {
                           RBETA = double( BETA );
                           BETA = CMPLX( RBETA, RZERO );
                        }
                        NULL = N <= 0;
                        if (CONJ) NULL = NULL || ( ( K <= 0 || ALPHA == ZERO ) && RBETA == RONE );

                        // Generate the matrix C.

                        cmake(SNAME( 2: 3 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1;

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO;
                        TRANSS = TRANS;
                        NS = N;
                        KS = K;
                        ALS = ALPHA;
                        for (I = 1; I <= LAA; I++) { // 10
                           AS[I] = AA( I );
                        } // 10
                        LDAS = LDA;
                        for (I = 1; I <= LBB; I++) { // 20
                           BS[I] = BB( I );
                        } // 20
                        LDBS = LDB;
                        if ( CONJ ) {
                           RBETS = RBETA;
                        } else {
                           BETS = BETA;
                        }
                        for (I = 1; I <= LCC; I++) { // 30
                           CS[I] = CC( I );
                        } // 30
                        LDCS = LDC;

                        // Call the subroutine.

                        if ( CONJ ) {
                           if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC;
                           if (REWI) REWIND NTRA;
                           cher2k(UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, RBETA, CC, LDC );
                        } else {
                           if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC;
                           if (REWI) REWIND NTRA;
                           csyr2k(UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9992 );
                           FATAL = true;
                           GO TO 150;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = UPLOS == UPLO;
                        ISAME[2] = TRANSS == TRANS;
                        ISAME[3] = NS == N;
                        ISAME[4] = KS == K;
                        ISAME[5] = ALS == ALPHA;
                        ISAME[6] = LCE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        ISAME[8] = LCE( BS, BB, LBB );
                        ISAME[9] = LDBS == LDB;
                        if ( CONJ ) {
                           ISAME[10] = RBETS == RBETA;
                        } else {
                           ISAME[10] = BETS == BETA;
                        }
                        if ( NULL ) {
                           ISAME[11] = LCE( CS, CC, LCC );
                        } else {
                           ISAME[11] = LCERES( 'HE', UPLO, N, N, CS, CC, LDC );
                        }
                        ISAME[12] = LDCS == LDC;

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME && ISAME( I );
                           if( !ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I;
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
                                    W[I] = ALPHA*AB( ( J - 1 )*2* NMAX + K + I );
                                    if ( CONJ ) {
                                       W[K + I] = CONJG( ALPHA )* AB( ( J - 1 )*2* NMAX + I );
                                    } else {
                                       W[K + I] = ALPHA* AB( ( J - 1 )*2* NMAX + I );
                                    }
                                 } // 50
                                 cmmch(TRANST, 'N', LJ, 1, 2*K, ONE, AB( JJAB ), 2*NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              } else {
                                 for (I = 1; I <= K; I++) { // 60
                                    if ( CONJ ) {
                                       W[I] = ALPHA*CONJG( AB( ( K + I - 1 )*NMAX + J ) )                                        W( K + I ) = CONJG( ALPHA* AB( ( I - 1 )*NMAX + J ) );
                                    } else {
                                       W[I] = ALPHA*AB( ( K + I - 1 )* NMAX + J )                                        W( K + I ) = ALPHA* AB( ( I - 1 )*NMAX + J );
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
                              ERRMAX = max( ERRMAX, ERR );
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
         WRITE( NOUT, FMT = 9999 )SNAME, NC;
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX;
      }
      GO TO 160;

      } // 140
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      if ( CONJ ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC;
      } else {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC;
      }

      } // 160
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), '(${.f4_1},${.f4_1}), A,${.i3}, B,${.i3},${.f4_1}, C,${.i3})           .' );
 9993 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), '(${.f4_1},${.f4_1}), A,${.i3}, B,${.i3},(${.f4_1},${.f4_1}), C,${.i3})    .' );
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void cchke(ISNUM,srnamc.SRNAMT, NOUT ) {

// Tests the error exits from the Level 3 Blas.
// Requires a special version of the error-handling routine XERBLA.
// A, B and C should not need to be defined.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

// 3-19-92:  Initialize ALPHA, BETA, RALPHA, and RBETA  (eca)
// 3-19-92:  Fix argument 12 in calls to CSYMM and CHEMM
             // with INFOT = 9  (eca)
      int                ISNUM, NOUT;
      String            srnamc.SRNAMT;
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Parameters ..
      double               ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      // .. Local Scalars ..
      Complex            ALPHA, BETA;
      double               RALPHA, RBETA;
      // .. Local Arrays ..
      Complex            A( 2, 1 ), B( 2, 1 ), C( 2, 1 );
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHEMM, CHER2K, CHERK, CHKXER, CSYMM, CSYR2K, CSYRK, CTRMM, CTRSM
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR

      // OK is set to false by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      OK = true;
      // LERR is set to true by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      LERR = false;

      // Initialize ALPHA, BETA, RALPHA, and RBETA.

      ALPHA = CMPLX( ONE, -ONE );
      BETA = CMPLX( TWO, -TWO );
      RALPHA = ONE;
      RBETA = TWO;

      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90 )ISNUM;
   10 INFOT = 1;
      cgemm('/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 1;
      cgemm('/', 'C', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 1;
      cgemm('/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemm('N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemm('C', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemm('T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('N', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('C', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('C', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('C', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('T', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemm('T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('N', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('C', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('C', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('C', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('T', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemm('T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('N', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('C', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('C', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('C', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('T', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemm('T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('N', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('C', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('C', 'C', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('C', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('T', 'C', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemm('T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('C', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('N', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('C', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('T', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('C', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemm('T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('N', 'C', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('C', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('C', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('C', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('T', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemm('T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   20 INFOT = 1;
      chemm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      chemm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      chemm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      chemm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      chemm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      chemm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      chemm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      chemm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      chemm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      chemm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      chemm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      chemm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      chemm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      chemm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      chemm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      chemm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      chemm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      chemm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      chemm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      chemm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      chemm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      chemm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   30 INFOT = 1;
      csymm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      csymm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csymm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csymm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csymm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csymm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csymm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csymm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csymm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csymm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csymm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csymm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csymm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csymm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csymm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csymm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csymm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csymm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csymm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csymm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csymm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csymm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   40 INFOT = 1;
      ctrmm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrmm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctrmm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrmm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrmm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrmm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'U', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'L', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrmm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrmm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   50 INFOT = 1;
      ctrsm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctrsm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctrsm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctrsm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctrsm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctrsm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'U', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'L', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctrsm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctrsm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   60 INFOT = 1;
      cherk('/', 'N', 0, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cherk('U', 'T', 0, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cherk('U', 'N', -1, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cherk('U', 'C', -1, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cherk('L', 'N', -1, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cherk('L', 'C', -1, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cherk('U', 'N', 0, -1, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cherk('U', 'C', 0, -1, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cherk('L', 'N', 0, -1, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cherk('L', 'C', 0, -1, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cherk('U', 'N', 2, 0, RALPHA, A, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cherk('U', 'C', 0, 2, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cherk('L', 'N', 2, 0, RALPHA, A, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cherk('L', 'C', 0, 2, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cherk('U', 'N', 2, 0, RALPHA, A, 2, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cherk('U', 'C', 2, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cherk('L', 'N', 2, 0, RALPHA, A, 2, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cherk('L', 'C', 2, 0, RALPHA, A, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   70 INFOT = 1;
      csyrk('/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      csyrk('U', 'C', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyrk('U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyrk('U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyrk('L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyrk('L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyrk('U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyrk('U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyrk('L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyrk('L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyrk('U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyrk('U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyrk('L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyrk('L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      csyrk('U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      csyrk('U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      csyrk('L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      csyrk('L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   80 INFOT = 1;
      cher2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cher2k('U', 'T', 0, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cher2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cher2k('U', 'C', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cher2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cher2k('L', 'C', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cher2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cher2k('U', 'C', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cher2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cher2k('L', 'C', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cher2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cher2k('U', 'C', 0, 2, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cher2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cher2k('L', 'C', 0, 2, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cher2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cher2k('U', 'C', 0, 2, ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cher2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, RBETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cher2k('L', 'C', 0, 2, ALPHA, A, 2, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cher2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cher2k('U', 'C', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cher2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cher2k('L', 'C', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 100;
   90 INFOT = 1;
      csyr2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      csyr2k('U', 'C', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyr2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyr2k('U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyr2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      csyr2k('L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyr2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyr2k('U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyr2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      csyr2k('L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyr2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyr2k('U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyr2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      csyr2k('L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csyr2k('U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      csyr2k('L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csyr2k('U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      csyr2k('L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );

  100 IF( OK )THEN;
         WRITE( NOUT, FMT = 9999 srnamc.SRNAMT;
      } else {
         WRITE( NOUT, FMT = 9998 srnamc.SRNAMT;
      }
      return;

 9999 FORMAT( ' ${.a6} PASSED THE TESTS OF ERROR-EXITS' );
 9998 FORMAT( ' ******* ${.a6} FAILED THE TESTS OF ERROR-EXITS *******' );
      }
      void cmake(final int TYPE, final int UPLO, final int DIAG, final int M, final int N, final int A, final int NMAX, final int AA, final int LDA, final int RESET, final int TRANSL,) {

// Generates values for an M by N matrix A.
// Stores the values in the array AA in the data structure required
// by the routine, with unwanted elements set to rogue value.

// TYPE is 'GE', 'HE', 'SY' or 'TR'.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e10, 1.0e10 ) ;
      double               RZERO;
      const              RZERO = 0.0 ;
      double               RROGUE;
      const              RROGUE = -1.0e10 ;
      // .. Scalar Arguments ..
      Complex            TRANSL;
      int                LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex            A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J, JJ;
      bool               GEN, HER, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      //- COMPLEX            CBEG;
      // EXTERNAL CBEG
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, REAL

      GEN = TYPE == 'GE';
      HER = TYPE == 'HE';
      SYM = TYPE == 'SY';
      TRI = TYPE == 'TR';
      UPPER = ( HER || SYM || TRI ) && UPLO == 'U';
      LOWER = ( HER || SYM || TRI ) && UPLO == 'L';
      UNIT = TRI && DIAG == 'U';

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN || ( UPPER && I <= J ) || ( LOWER && I >= J ) ) {
               A[I][J] = CBEG( RESET ) + TRANSL;
               if ( I != J ) {
                  // Set some elements to zero
                  if (N > 3 && J == N/2) A( I, J ) = ZERO;
                  if ( HER ) {
                     A[J][I] = CONJG( A( I, J ) );
                  } else if ( SYM ) {
                     A[J][I] = A( I, J );
                  } else if ( TRI ) {
                     A[J][I] = ZERO;
                  }
               }
            }
         } // 10
         if (HER) A( J, J ) = CMPLX( double( A( J, J ) ), RZERO );
         if[TRI ) A( J][J] = A( J, J ) + ONE;
         IF[UNIT ) A( J][J] = ONE;
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA[I + ( J - 1 )*LDA] = A( I, J );
            } // 30
            for (I = M + 1; I <= LDA; I++) { // 40
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 40
         } // 50
      } else if ( TYPE == 'HE' || TYPE == 'SY' || TYPE == 'TR' ) {
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
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 60
            for (I = IBEG; I <= IEND; I++) { // 70
               AA[I + ( J - 1 )*LDA] = A( I, J );
            } // 70
            for (I = IEND + 1; I <= LDA; I++) { // 80
               AA[I + ( J - 1 )*LDA] = ROGUE;
            } // 80
            if ( HER ) {
               JJ = J + ( J - 1 )*LDA;
               AA[JJ] = CMPLX( double( AA( JJ ) ), RROGUE );
            }
         } // 90
      }
      return;
      }
      void cmmch(final int TRANSA, final int TRANSB, final int M, final int N, final int KK, final int ALPHA, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int BETA, final Matrix<double> C_, final int LDC, final int CT, final int G, final Matrix<double> CC_, final int LDCC, final int EPS, final int ERR, final int FATAL, final int NOUT, final int MV,) {
  final A = A_.dim();
  final B = B_.dim();
  final C = C_.dim();
  final CC = CC_.dim();

// Checks the results of the computational tests.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      double               RZERO, RONE;
      const              RZERO = 0.0, RONE = 1.0 ;
      // .. Scalar Arguments ..
      Complex            ALPHA, BETA;
      double               EPS, ERR;
      int                KK, LDA, LDB, LDC, LDCC, M, N, NOUT;
      bool               FATAL, MV;
      String             TRANSA, TRANSB;
      // .. Array Arguments ..
      Complex            A( LDA, * ), B( LDB, * ), C( LDC, * ), CC( LDCC, * ), CT( * );
      double               G( * );
      // .. Local Scalars ..
      Complex            CL;
      double               ERRI;
      int                I, J, K;
      bool               CTRANA, CTRANB, TRANA, TRANB;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
      // .. Statement Functions ..
      double               ABS1;
      // .. Statement Function definitions ..
      ABS1[CL] = ( double( CL ) ).abs() + ( AIMAG( CL ) ).abs();

      TRANA = TRANSA == 'T' || TRANSA == 'C';
      TRANB = TRANSB == 'T' || TRANSB == 'C';
      CTRANA = TRANSA == 'C';
      CTRANB = TRANSB == 'C';

      // Compute expected result, one column at a time, in CT using data
      // in A, B and C.
      // Compute gauges in G.

      for (J = 1; J <= N; J++) { // 220

         for (I = 1; I <= M; I++) { // 10
            CT[I] = ZERO;
            G[I] = RZERO;
         } // 10
         if ( !TRANA && !TRANB ) {
            for (K = 1; K <= KK; K++) { // 30
               for (I = 1; I <= M; I++) { // 20
                  CT[I] = CT( I ) + A( I, K )*B( K, J );
                  G[I] = G( I ) + ABS1( A( I, K ) )*ABS1( B( K, J ) );
               } // 20
            } // 30
         } else if ( TRANA && !TRANB ) {
            if ( CTRANA ) {
               for (K = 1; K <= KK; K++) { // 50
                  for (I = 1; I <= M; I++) { // 40
                     CT[I] = CT( I ) + CONJG( A( K, I ) )*B( K, J );
                     G[I] = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) );
                  } // 40
               } // 50
            } else {
               for (K = 1; K <= KK; K++) { // 70
                  for (I = 1; I <= M; I++) { // 60
                     CT[I] = CT( I ) + A( K, I )*B( K, J );
                     G[I] = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) );
                  } // 60
               } // 70
            }
         } else if ( !TRANA && TRANB ) {
            if ( CTRANB ) {
               for (K = 1; K <= KK; K++) { // 90
                  for (I = 1; I <= M; I++) { // 80
                     CT[I] = CT( I ) + A( I, K )*CONJG( B( J, K ) );
                     G[I] = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) );
                  } // 80
               } // 90
            } else {
               for (K = 1; K <= KK; K++) { // 110
                  for (I = 1; I <= M; I++) { // 100
                     CT[I] = CT( I ) + A( I, K )*B( J, K );
                     G[I] = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) );
                  } // 100
               } // 110
            }
         } else if ( TRANA && TRANB ) {
            if ( CTRANA ) {
               if ( CTRANB ) {
                  for (K = 1; K <= KK; K++) { // 130
                     for (I = 1; I <= M; I++) { // 120
                        CT[I] = CT( I ) + CONJG( A( K, I ) )* CONJG( B( J, K ) )                         G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 120
                  } // 130
               } else {
                  for (K = 1; K <= KK; K++) { // 150
                     for (I = 1; I <= M; I++) { // 140
                        CT[I] = CT( I ) + CONJG( A( K, I ) )*B( J, K );
                        G[I] = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 140
                  } // 150
               }
            } else {
               if ( CTRANB ) {
                  for (K = 1; K <= KK; K++) { // 170
                     for (I = 1; I <= M; I++) { // 160
                        CT[I] = CT( I ) + A( K, I )*CONJG( B( J, K ) );
                        G[I] = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 160
                  } // 170
               } else {
                  for (K = 1; K <= KK; K++) { // 190
                     for (I = 1; I <= M; I++) { // 180
                        CT[I] = CT( I ) + A( K, I )*B( J, K );
                        G[I] = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) );
                     } // 180
                  } // 190
               }
            }
         }
         for (I = 1; I <= M; I++) { // 200
            CT[I] = ALPHA*CT( I ) + BETA*C( I, J );
            G[I] = ABS1( ALPHA )*G( I ) + ABS1( BETA )*ABS1( C( I, J ) );
         } // 200

         // Compute the error ratio for this result.

         ERR = ZERO;
         for (I = 1; I <= M; I++) { // 210
            ERRI = ABS1( CT( I ) - CC( I, J ) )/EPS;
            if( G( I ) != RZERO ) ERRI = ERRI/G( I );
            ERR = max( ERR, ERRI );
            if( ERR*sqrt( EPS ) >= RONE ) GO TO 230;
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

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT' );
 9998 FORMAT( 1X, I7, 2( '  (${.g15_6},${.g15_6})' ) );
 9997 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
      }
      bool lce(final int RI, final int RJ, final int LR,) {

// Tests if two arrays are identical.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      int                LR;
      // .. Array Arguments ..
      Complex            RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;

      for (I = 1; I <= LR; I++) { // 10
         if( RI( I ) != RJ( I ) ) GO TO 20;
      } // 10
      LCE = true;
      GO TO 30;
      } // 20
      LCE = false;
   30 return;
      }
      bool lceres(final int TYPE, final int UPLO, final int M, final int N, final int AA, final int AS, final int LDA,) {

// Tests if selected elements in two arrays are equal.

// TYPE is 'GE' or 'HE' or 'SY'.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      Complex            AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;

      UPPER = UPLO == 'U';
      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 10
         } // 20
      } else if ( TYPE == 'HE' || TYPE == 'SY' ) {
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

      LCERES = true;
      GO TO 80;
      } // 70
      LCERES = false;
   80 return;
      }
      Complex cbeg(final int RESET,) {

// Generates complex numbers as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ;
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX

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
      }
      double sdiff(final int X, final int Y,) {

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      double               X, Y;

      SDIFF = X - Y;
      return;
      }
      void chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK ) {

// Tests whether XERBLA has detected an error when it should.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      int                INFOT, NOUT;
      bool               LERR, OK;
      String            srnamc.SRNAMT;

      if ( !LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT,srnamc.SRNAMT;
         OK = false;
      }
      LERR = false;
      return;

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${.i2} NOT DETECTED BY ${.a6} *****' );
      }
      void xerbla(final int SRNAME, final Box<int> INFO,) {

// This is a special version of XERBLA to be used only as part of
// the test program for testing error exits from the Level 3 BLAS
// routines.

// XERBLA  is an error handler for the Level 3 BLAS routines.

// It is called by the Level 3 BLAS routines if an input parameter is
// invalid.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      int                INFO;
      String             SRNAME;
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUT, OK, LERR
      // COMMON /SRNAMC(srnamc.SRNAMT

      LERR = true;
      if ( INFO != INFOT ) {
         if ( INFOT != 0 ) {
            WRITE( NOUT, FMT = 9999 )INFO, INFOT;
         } else {
            WRITE( NOUT, FMT = 9997 )INFO;
         }
         OK = false;
      }
      if ( SRNAME !=srnamc.SRNAMT ) {
         WRITE( NOUT, FMT = 9998 )SRNAME,srnamc.SRNAMT;
         OK = false;
      }
      return;

 9999 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ${.i6} INSTEAD OF ${.i2} *******' );
 9998 FORMAT( ' ******* XERBLA WAS CALLED WITH SRNAME = ${.a6} INSTEAD OF ${.a6} *******' );
 9997 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ${.i6} *******' );
      }