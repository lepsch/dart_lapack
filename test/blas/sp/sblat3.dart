      void main() {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 6 ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      double               AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBEMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX ), W( 2*NMAX );
      int                IDIM( NIDMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      //- REAL               SDIFF;
      //- bool               LSE;
      // EXTERNAL SDIFF, LSE
      // .. External Subroutines ..
      // EXTERNAL SCHK1, SCHK2, SCHK3, SCHK4, SCHK5, SCHKE, SMMCH
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
      const SNAMES = ['SGEMM ', 'SSYMM ', 'STRMM ', 'STRSM ', 'SSYRK ', 'SSYR2K'];
      // .. Executable Statements ..

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

      // } // 60
      CLOSE ( NIN );

      // Compute EPS (the machine precision).

      EPS = EPSILON(ZERO);
      WRITE( NOUT, FMT = 9998 )EPS;

      // Check the reliability of SMMCH using exact data.

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
      // CC holds the exact result. On exit from SMMCH CT holds
      // the result computed by SMMCH.
      TRANSA = 'N';
      TRANSB = 'N';
      smmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( CC, CT, N );
      if ( !SAME || ERR != ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }
      TRANSB = 'T';
      smmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( CC, CT, N );
      if ( !SAME || ERR != ZERO ) {
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
      TRANSA = 'T';
      TRANSB = 'N';
      smmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( CC, CT, N );
      if ( !SAME || ERR != ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR;
         STOP;
      }
      TRANSB = 'T';
      smmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( CC, CT, N );
      if ( !SAME || ERR != ZERO ) {
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
               schke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * );
            }
            // Test computations.
            INFOT = 0;
            OK = true;
            FATAL = false;
            GO TO ( 140, 150, 160, 160, 170, 180 )ISNUM;
            // Test SGEMM, 01.
  // 140       CALL SCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test SSYMM, 02.
  // 150       CALL SCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test STRMM, 03, STRSM, 04.
  // 160       CALL SCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C );
            GO TO 190;
            // Test SSYRK, 05.
  // 170       CALL SCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G );
            GO TO 190;
            // Test SSYR2K, 06.
  // 180       CALL SCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W );
            GO TO 190;

  // 190       IF( FATAL && SFATAL ) GO TO 210;
         }
      } // 200
      WRITE( NOUT, FMT = 9986 );
      GO TO 230;

      // } // 210
      WRITE( NOUT, FMT = 9985 );
      GO TO 230;

      // } // 220
      WRITE( NOUT, FMT = 9991 );

      // } // 230
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT );
      STOP;

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LESS THAN${.f8_2}');
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 );
 9997 FORMAT( ' NUMBER OF VALUES OF ${} IS LESS THAN 1 OR GREATER THAN ${.i2}');
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${.i2}');
 9995 FORMAT( ' TESTS OF THE REAL             LEVEL 3 BLAS', //' THE FOLLOWING PARAMETER VALUES WILL BE USED:' );
 9994 FORMAT( '   FOR N              ', 9I6 );
 9993 FORMAT( '   FOR ALPHA          ', 7F6.1 );
 9992 FORMAT( '   FOR BETA           ', 7F6.1 );
 9991 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******' );
 9990 FORMAT( ' SUBPROGRAM NAME ${.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******' );
 9989 FORMAT( ' ERROR IN SMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n SMMCH WAS CALLED WITH TRANSA = ${.a1} AND TRANSB = ${.a1}\n AND RETURNED SAME = ${.l1} AND ERR = ${.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******' );
 9988 FORMAT( A6, L2 );
 9987 FORMAT(' ${.a6} WAS NOT TESTED' );
 9986 FORMAT('\n END OF TESTS' );
 9985 FORMAT('\n ******* FATAL ERROR - TESTS ABANDONED *******' );
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' );
      }

      void schk1(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G ) {

// Tests SGEMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double               ALPHA, ALS, BETA, BLS, ERR, ERRMAX;
      int                I, IA, IB, ICA, ICB, IK, IM, IN, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MA, MB, MS, N, NA, NARGS, NB, NC, NS;
      bool               NULL, RESET, SAME, TRANA, TRANB;
      String             TRANAS, TRANBS, TRANSA, TRANSB;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SMAKE, SMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICH = 'NTC';
      // .. Executable Statements ..

      NARGS = 13;
      NC = 0;
      RESET = true;
      ERRMAX = ZERO;

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

                  smake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

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

                     smake('GE', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO );

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA );

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB );

                           // Generate the matrix C.

                           smake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                           sgemm(TRANSA, TRANSB, M, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

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
                           ISAME[7] = LSE( AS, AA, LAA );
                           ISAME[8] = LDAS == LDA;
                           ISAME[9] = LSE( BS, BB, LBB );
                           ISAME[10] = LDBS == LDB;
                           ISAME[11] = BLS == BETA;
                           if ( NULL ) {
                              ISAME[12] = LSE( CS, CC, LCC );
                           } else {
                              ISAME[12] = LSERES( 'GE', ' ', M, N, CS, CC, LDC );
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

                              smmch(TRANSA, TRANSB, M, N, K, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
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

      // } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC;

      // } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(''${.a1}'',''${.a1}'',${i3(3, ',')}${.f4_1}, A,${.i3}, B,${.i3},${.f4_1}, C,${.i3}).' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }

      void schk2(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G ) {

// Tests SSYMM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double               ALPHA, ALS, BETA, BLS, ERR, ERRMAX;
      int                I, IA, IB, ICS, ICU, IM, IN, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String             SIDE, SIDES, UPLO, UPLOS;
      String             ICHS, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMMCH, SSYMM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHS = 'LR', ICHU = 'UL';
      // .. Executable Statements ..

      NARGS = 12;
      NC = 0;
      RESET = true;
      ERRMAX = ZERO;

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

            smake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

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

                  // Generate the symmetric matrix A.

                  smake('SY', UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );

                        // Generate the matrix C.

                        smake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                        ssymm(SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

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
                        ISAME[6] = LSE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        ISAME[8] = LSE( BS, BB, LBB );
                        ISAME[9] = LDBS == LDB;
                        ISAME[10] = BLS == BETA;
                        if ( NULL ) {
                           ISAME[11] = LSE( CS, CC, LCC );
                        } else {
                           ISAME[11] = LSERES( 'GE', ' ', M, N, CS, CC, LDC );
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
                              smmch('N', 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
                           } else {
                              smmch('N', 'N', M, N, N, ALPHA, B, NMAX, A, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, true );
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

      // } // 110
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC;

      // } // 120
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), F4.1, ', A,${.i3}, B,${.i3},${.f4_1}, C,${.i3})    .' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }


      void schk3(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, A, AA, AS, B, BB, BS, CT, G, C ) {

// Tests STRMM and STRSM.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double               ALPHA, ALS, ERR, ERRMAX;
      int                I, IA, ICD, ICS, ICT, ICU, IM, IN, J, LAA, LBB, LDA, LDAS, LDB, LDBS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String             DIAG, DIAGS, SIDE, SIDES, TRANAS, TRANSA, UPLO, UPLOS;
      String             ICHD, ICHS, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMMCH, STRMM, STRSM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN', ICHS = 'LR';
      // .. Executable Statements ..

      NARGS = 11;
      NC = 0;
      RESET = true;
      ERRMAX = ZERO;
      // Set up zero matrix for SMMCH.
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

                           smake('TR', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                           // Generate the matrix B.

                           smake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

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
                              strmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           } else if ( SNAME( 4: 5 ) == 'SM' ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;
                              if (REWI) REWIND NTRA;
                              strsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
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
                           ISAME[8] = LSE( AS, AA, LAA );
                           ISAME[9] = LDAS == LDA;
                           if ( NULL ) {
                              ISAME[10] = LSE( BS, BB, LBB );
                           } else {
                              ISAME[10] = LSERES( 'GE', ' ', M, N, BS, BB, LDB );
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
                                    smmch(TRANSA, 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
                                 } else {
                                    smmch('N', TRANSA, M, N, N, ALPHA, B, NMAX, A, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, true );
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
                                    smmch(TRANSA, 'N', M, N, M, ONE, A, NMAX, C, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
                                 } else {
                                    smmch('N', TRANSA, M, N, N, ONE, C, NMAX, A, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, false );
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

      // } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;

      // } // 160
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT(' ${.i6}: ${.a6}(', 4( '\'${.a1}\',' ), 2( I3, ',' ), F4.1, ', A,${.i3}, B,${.i3})        .' );
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void schk4(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G ) {

// Tests SSYRK.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double               ALPHA, ALS, BETA, BETS, ERR, ERRMAX;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, KS, LAA, LCC, LDA, LDAS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, UPLO, UPLOS;
      String             ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMMCH, SSYRK
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHT = 'NTC', ICHU = 'UL';
      // .. Executable Statements ..

      NARGS = 10;
      NC = 0;
      RESET = true;
      ERRMAX = ZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN );
         // Set LDC to 1 more than minimum value if room.
         LDC = N;
         if (LDC < NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC > NMAX) GO TO 100;
         LCC = LDC*N;
         NULL = N <= 0;

         for (IK = 1; IK <= NIDIM; IK++) { // 90
            K = IDIM( IK );

            for (ICT = 1; ICT <= 3; ICT++) { // 80
               TRANS = ICHT( ICT: ICT );
               TRAN = TRANS == 'T' || TRANS == 'C';
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

               smake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB );

                        // Generate the matrix C.

                        smake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                        BETS = BETA;
                        for (I = 1; I <= LCC; I++) { // 20
                           CS[I] = CC( I );
                        } // 20
                        LDCS = LDC;

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        ssyrk(UPLO, TRANS, N, K, ALPHA, AA, LDA, BETA, CC, LDC );

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9993 );
                           FATAL = true;
                           GO TO 120;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = UPLOS == UPLO;
                        ISAME[2] = TRANSS == TRANS;
                        ISAME[3] = NS == N;
                        ISAME[4] = KS == K;
                        ISAME[5] = ALS == ALPHA;
                        ISAME[6] = LSE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        ISAME[8] = BETS == BETA;
                        if ( NULL ) {
                           ISAME[9] = LSE( CS, CC, LCC );
                        } else {
                           ISAME[9] = LSERES( 'SY', UPLO, N, N, CS, CC, LDC );
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
                                 smmch('T', 'N', LJ, 1, K, ALPHA, A( 1, JJ ), NMAX, A( 1, J ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              } else {
                                 smmch('N', 'T', LJ, 1, K, ALPHA, A( JJ, 1 ), NMAX, A( J, 1 ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
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

      // } // 110
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      // } // 120
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC;

      // } // 130
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), F4.1, ', A,${.i3},${.f4_1}, C,${.i3})           .' );
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void schk5(SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W ) {

// Tests SSYR2K.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // .. Scalar Arguments ..
      double               EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double               AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX ), W( 2*NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double               ALPHA, ALS, BETA, BETS, ERR, ERRMAX;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, UPLO, UPLOS;
      String             ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      //- bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMMCH, SSYR2K
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      const ICHT = 'NTC', ICHU = 'UL';
      // .. Executable Statements ..

      NARGS = 12;
      NC = 0;
      RESET = true;
      ERRMAX = ZERO;

      for (IN = 1; IN <= NIDIM; IN++) { // 130
         N = IDIM( IN );
         // Set LDC to 1 more than minimum value if room.
         LDC = N;
         if (LDC < NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC > NMAX) GO TO 130;
         LCC = LDC*N;
         NULL = N <= 0;

         for (IK = 1; IK <= NIDIM; IK++) { // 120
            K = IDIM( IK );

            for (ICT = 1; ICT <= 3; ICT++) { // 110
               TRANS = ICHT( ICT: ICT );
               TRAN = TRANS == 'T' || TRANS == 'C';
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
                  smake('GE', ' ', ' ', MA, NA, AB, 2*NMAX, AA, LDA, RESET, ZERO );
               } else {
                  smake('GE', ' ', ' ', MA, NA, AB, NMAX, AA, LDA, RESET, ZERO );
               }

               // Generate the matrix B.

               LDB = LDA;
               LBB = LAA;
               if ( TRAN ) {
                  smake('GE', ' ', ' ', MA, NA, AB( K + 1 ), 2*NMAX, BB, LDB, RESET, ZERO );
               } else {
                  smake('GE', ' ', ' ', MA, NA, AB( K*NMAX + 1 ), NMAX, BB, LDB, RESET, ZERO );
               }

               for (ICU = 1; ICU <= 2; ICU++) { // 100
                  UPLO = ICHU( ICU: ICU );
                  UPPER = UPLO == 'U';

                  for (IA = 1; IA <= NALF; IA++) { // 90
                     ALPHA = ALF( IA );

                     for (IB = 1; IB <= NBET; IB++) { // 80
                        BETA = BET( IB );

                        // Generate the matrix C.

                        smake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

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
                        BETS = BETA;
                        for (I = 1; I <= LCC; I++) { // 30
                           CS[I] = CC( I );
                        } // 30
                        LDCS = LDC;

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        ssyr2k(UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                        // Check if error-exit was taken incorrectly.

                        if ( !OK ) {
                           WRITE( NOUT, FMT = 9993 );
                           FATAL = true;
                           GO TO 150;
                        }

                        // See what data changed inside subroutines.

                        ISAME[1] = UPLOS == UPLO;
                        ISAME[2] = TRANSS == TRANS;
                        ISAME[3] = NS == N;
                        ISAME[4] = KS == K;
                        ISAME[5] = ALS == ALPHA;
                        ISAME[6] = LSE( AS, AA, LAA );
                        ISAME[7] = LDAS == LDA;
                        ISAME[8] = LSE( BS, BB, LBB );
                        ISAME[9] = LDBS == LDB;
                        ISAME[10] = BETS == BETA;
                        if ( NULL ) {
                           ISAME[11] = LSE( CS, CC, LCC );
                        } else {
                           ISAME[11] = LSERES( 'SY', UPLO, N, N, CS, CC, LDC );
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
                                    W[I] = AB( ( J - 1 )*2*NMAX + K + I )                                     W( K + I ) = AB( ( J - 1 )*2*NMAX + I );
                                 } // 50
                                 smmch('T', 'N', LJ, 1, 2*K, ALPHA, AB( JJAB ), 2*NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
                              } else {
                                 for (I = 1; I <= K; I++) { // 60
                                    W[I] = AB( ( K + I - 1 )*NMAX + J )                                     W( K + I ) = AB( ( I - 1 )*NMAX + J );
                                 } // 60
                                 smmch('N', 'N', LJ, 1, 2*K, ALPHA, AB( JJ ), NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, true );
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

      // } // 140
      if (N > 1) WRITE( NOUT, FMT = 9995 )J;

      // } // 150
      WRITE( NOUT, FMT = 9996 )SNAME;
      WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC;

      // } // 160
      return;

 9999 FORMAT( ' ${.a6} PASSED THE COMPUTATIONAL TESTS (${.i6} CALLS)' );
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ${.i2} WAS CHANGED INCORRECTLY *******' );
 9997 FORMAT( ' ${.a6} COMPLETED THE COMPUTATIONAL TESTS (${.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${.f8_2} - SUSPECT *******' );
 9996 FORMAT( ' ******* ${.a6} FAILED ON CALL NUMBER:' );
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
 9994 FORMAT(' ${.i6}: ${.a6}(', 2( '\'${.a1}\',' ), 2( I3, ',' ), F4.1, ', A,${.i3}, B,${.i3},${.f4_1}, C,${.i3})    .' );
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******' );
      }
      void schke(ISNUM,srnamc.SRNAMT, NOUT ) {

// Tests the error exits from the Level 3 Blas.
// Requires a special version of the error-handling routine XERBLA.
// A, B and C should not need to be defined.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

// 3-19-92:  Initialize ALPHA and BETA  (eca)
// 3-19-92:  Fix argument 12 in calls to SSYMM with INFOT = 9  (eca)
      int                ISNUM, NOUT;
      String            srnamc.SRNAMT;
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Parameters ..
      double               ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      // .. Local Scalars ..
      double               ALPHA, BETA;
      // .. Local Arrays ..
      double               A( 2, 1 ), B( 2, 1 ), C( 2, 1 );
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGEMM, SSYMM, SSYR2K, SSYRK, STRMM, STRSM
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      // OK is set to false by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      OK = true;
      // LERR is set to true by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      LERR = false;

      // Initialize ALPHA and BETA.

      ALPHA = ONE;
      BETA = TWO;

      GO TO ( 10, 20, 30, 40, 50, 60 )ISNUM;
   10 INFOT = 1;
      sgemm('/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 1;
      sgemm('/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemm('N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemm('T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemm('N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemm('N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemm('T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemm('T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemm('N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemm('N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemm('T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemm('T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemm('N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemm('N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemm('T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemm('T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemm('N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemm('N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemm('T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemm('T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemm('N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemm('T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemm('N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemm('T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemm('N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemm('N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemm('T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemm('T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70;
   20 INFOT = 1;
      ssymm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ssymm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssymm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssymm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssymm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssymm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssymm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssymm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssymm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssymm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssymm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssymm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssymm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssymm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssymm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssymm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssymm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssymm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssymm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssymm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssymm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssymm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70;
   30 INFOT = 1;
      strmm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      strmm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      strmm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      strmm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strmm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strmm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strmm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strmm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70;
   40 INFOT = 1;
      strsm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      strsm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      strsm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      strsm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strsm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6;
      strsm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      strsm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11;
      strsm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70;
   50 INFOT = 1;
      ssyrk('/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ssyrk('U', '/', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyrk('U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyrk('U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyrk('L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyrk('L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyrk('U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyrk('U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyrk('L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyrk('L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyrk('U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyrk('U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyrk('L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyrk('L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ssyrk('U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ssyrk('U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ssyrk('L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ssyrk('L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70;
   60 INFOT = 1;
      ssyr2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ssyr2k('U', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyr2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyr2k('U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyr2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssyr2k('L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyr2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyr2k('U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyr2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssyr2k('L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyr2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyr2k('U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyr2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ssyr2k('L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssyr2k('U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ssyr2k('L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssyr2k('U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12;
      ssyr2k('L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(srnamc.SRNAMT, INFOT, NOUT, LERR, OK );

      // } 70
      if( OK ) {
         WRITE( NOUT, FMT = 9999 srnamc.SRNAMT;
      } else {
         WRITE( NOUT, FMT = 9998 srnamc.SRNAMT;
      }
      return;

 9999 FORMAT( ' ${.a6} PASSED THE TESTS OF ERROR-EXITS' );
 9998 FORMAT( ' ******* ${.a6} FAILED THE TESTS OF ERROR-EXITS *******' );
      }
      void smake(TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, RESET, TRANSL ) {

// Generates values for an M by N matrix A.
// Stores the values in the array AA in the data structure required
// by the routine, with unwanted elements set to rogue value.

// TYPE is 'GE', 'SY' or 'TR'.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ROGUE;
      const              ROGUE = -1.0e10 ;
      // .. Scalar Arguments ..
      double               TRANSL;
      int                LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      double               A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      //- REAL               SBEG;
      // EXTERNAL SBEG
      // .. Executable Statements ..
      GEN = TYPE == 'GE';
      SYM = TYPE == 'SY';
      TRI = TYPE == 'TR';
      UPPER = ( SYM || TRI ) && UPLO == 'U';
      LOWER = ( SYM || TRI ) && UPLO == 'L';
      UNIT = TRI && DIAG == 'U';

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN || ( UPPER && I <= J ) || ( LOWER && I >= J ) ) {
               A[I][J] = SBEG( RESET ) + TRANSL;
               if ( I != J ) {
                  // Set some elements to zero
                  if (N > 3 && J == N/2) A( I, J ) = ZERO;
                  if ( SYM ) {
                     A[J][I] = A( I, J );
                  } else if ( TRI ) {
                     A[J][I] = ZERO;
                  }
               }
            }
         } // 10
         if (TRI) A( J, J ) = A( J, J ) + ONE;
         if(UNIT ) A[J][J] = ONE;
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
      } else if ( TYPE == 'SY' || TYPE == 'TR' ) {
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
         } // 90
      }
      return;
      }
      void smmch(TRANSA, TRANSB, M, N, KK, ALPHA, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, BETA, final Matrix<double> C, final int LDC, CT, G, final Matrix<double> CC, final int LDCC, EPS, ERR, FATAL, NOUT, MV ) {

// Checks the results of the computational tests.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // .. Scalar Arguments ..
      double               ALPHA, BETA, EPS, ERR;
      int                KK, LDA, LDB, LDC, LDCC, M, N, NOUT;
      bool               FATAL, MV;
      String             TRANSA, TRANSB;
      // .. Array Arguments ..
      double               A( LDA, * ), B( LDB, * ), C( LDC, * ), CC( LDCC, * ), CT( * ), G( * );
      // .. Local Scalars ..
      double               ERRI;
      int                I, J, K;
      bool               TRANA, TRANB;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // .. Executable Statements ..
      TRANA = TRANSA == 'T' || TRANSA == 'C';
      TRANB = TRANSB == 'T' || TRANSB == 'C';

      // Compute expected result, one column at a time, in CT using data
      // in A, B and C.
      // Compute gauges in G.

      for (J = 1; J <= N; J++) { // 120

         for (I = 1; I <= M; I++) { // 10
            CT[I] = ZERO;
            G[I] = ZERO;
         } // 10
         if ( !TRANA && !TRANB ) {
            for (K = 1; K <= KK; K++) { // 30
               for (I = 1; I <= M; I++) { // 20
                  CT[I] = CT( I ) + A( I, K )*B( K, J );
                  G[I] = G( I ) + ( A( I, K ) ).abs()*( B( K, J ) ).abs();
               } // 20
            } // 30
         } else if ( TRANA && !TRANB ) {
            for (K = 1; K <= KK; K++) { // 50
               for (I = 1; I <= M; I++) { // 40
                  CT[I] = CT( I ) + A( K, I )*B( K, J );
                  G[I] = G( I ) + ( A( K, I ) ).abs()*( B( K, J ) ).abs();
               } // 40
            } // 50
         } else if ( !TRANA && TRANB ) {
            for (K = 1; K <= KK; K++) { // 70
               for (I = 1; I <= M; I++) { // 60
                  CT[I] = CT( I ) + A( I, K )*B( J, K );
                  G[I] = G( I ) + ( A( I, K ) ).abs()*( B( J, K ) ).abs();
               } // 60
            } // 70
         } else if ( TRANA && TRANB ) {
            for (K = 1; K <= KK; K++) { // 90
               for (I = 1; I <= M; I++) { // 80
                  CT[I] = CT( I ) + A( K, I )*B( J, K );
                  G[I] = G( I ) + ( A( K, I ) ).abs()*( B( J, K ) ).abs();
               } // 80
            } // 90
         }
         for (I = 1; I <= M; I++) { // 100
            CT[I] = ALPHA*CT( I ) + BETA*C( I, J );
            G[I] = ( ALPHA ).abs()*G( I ) + ( BETA ).abs()*( C( I, J ) ).abs();
         } // 100

         // Compute the error ratio for this result.

         ERR = ZERO;
         for (I = 1; I <= M; I++) { // 110
            ERRI = ABS( CT( I ) - CC( I, J ) )/EPS;
            if( G( I ) != ZERO ) ERRI = ERRI/G( I );
            ERR = max( ERR, ERRI );
            if( ERR*sqrt( EPS ) >= ONE ) GO TO 130;
         } // 110

      } // 120

      // If the loop completes, all results are at least half accurate.
      GO TO 150;

      // Report fatal error.

  130 FATAL = true;
      WRITE( NOUT, FMT = 9999 );
      for (I = 1; I <= M; I++) { // 140
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, CT( I ), CC( I, J );
         } else {
            WRITE( NOUT, FMT = 9998 )I, CC( I, J ), CT( I );
         }
      } // 140
      if (N > 1) WRITE( NOUT, FMT = 9997 )J;

      // } // 150
      return;

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT' );
 9998 FORMAT( 1X, I7, 2G18.6 );
 9997 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ${.i3}');
      }
      bool lse(RI, RJ, LR ) {

// Tests if two arrays are identical.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      int                LR;
      // .. Array Arguments ..
      double               RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         if( RI( I ) != RJ( I ) ) GO TO 20;
      } // 10
      LSE = true;
      GO TO 30;
      // } // 20
      LSE = false;
   30 return;
      }
      bool lseres(TYPE, UPLO, M, N, AA, AS, LDA ) {

// Tests if selected elements in two arrays are equal.

// TYPE is 'GE' or 'SY'.

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
      double               AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO == 'U';
      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               if( AA( I, J ) != AS( I, J ) ) GO TO 70;
            } // 10
         } // 20
      } else if ( TYPE == 'SY' ) {
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

      LSERES = true;
      GO TO 80;
      // } // 70xx
      LSERES = false;
   80 return;
      }
      double sbeg(RESET ) {

// Generates random numbers uniformly distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, MI;
      // .. Save statement ..
      SAVE               I, IC, MI;
      // .. Executable Statements ..
    if ( RESET ) {
         // Initialize local variables.
         MI = 891;
         I = 7;
         IC = 0;
         RESET = false;
      }

      // The sequence of values of I is bounded between 1 and 999.
      // If initial I = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I = 4 or 8, the period will be 25.
      // If initial I = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I in 6.

      IC = IC + 1;
   10 I = I*MI;
      I = I - 1000*( I/1000 );
      if ( IC >= 5 ) {
         IC = 0;
         GO TO 10;
      }
      SBEG = ( I - 500 )/1001.0;
      return;
      }
      double sdiff(X, Y ) {

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.
      double               X, Y;
      // .. Executable Statements ..
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
      // .. Executable Statements ..
      if ( !LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT,srnamc.SRNAMT;
         OK = false;
      }
      LERR = false;
      return;

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${.i2} NOT DETECTED BY ${.a6} *****' );
      }
      void xerbla(SRNAME, Box<int> INFO ) {

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
      // .. Executable Statements ..
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
