void main() {
*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 6 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      int                NMAX;
      const              NMAX = 65 ;
      int                NIDMAX, NALMAX, NBEMAX;
      const              NIDMAX = 9, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      double             EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NOUT, NTRA;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR;
      String             TRANSA, TRANSB;
      String             SNAMET;
      String             SNAPS, SUMMRY;
      // .. Local Arrays ..
      double             AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBEMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX ), W( 2*NMAX );
      int                IDIM( NIDMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      double             DDIFF;
      bool               LDE;
      // EXTERNAL DDIFF, LDE
      // .. External Subroutines ..
      // EXTERNAL DCHK1, DCHK2, DCHK3, DCHK4, DCHK5, DCHKE, DMMCH
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
      DATA               SNAMES/'DGEMM ', 'DSYMM ', 'DTRMM ', 'DTRSM ', 'DSYRK ', 'DSYR2K'/
      // .. Executable Statements ..

      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SUMMRY
      READ( NIN, FMT = * )NOUT
      OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
      NOUTC = NOUT

      // Read name and unit number for snapshot output file and open file.

      READ( NIN, FMT = * )SNAPS
      READ( NIN, FMT = * )NTRA
      TRACE = NTRA.GE.0
      if ( TRACE ) {
         OPEN( NTRA, FILE = SNAPS, STATUS = 'UNKNOWN' )
      }
      // Read the flag that directs rewinding of the snapshot file.
      READ( NIN, FMT = * )REWI
      REWI = REWI.AND.TRACE
      // Read the flag that directs stopping on any failure.
      READ( NIN, FMT = * )SFATAL
      // Read the flag that indicates whether error exits are to be tested.
      READ( NIN, FMT = * )TSTERR
      // Read the threshold value of the test ratio
      READ( NIN, FMT = * )THRESH

      // Read and check the parameter values for the tests.

      // Values of N
      READ( NIN, FMT = * )NIDIM
      if ( NIDIM.LT.1.OR.NIDIM.GT.NIDMAX ) {
         WRITE( NOUT, FMT = 9997 )'N', NIDMAX
         GO TO 220
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM )
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ).LT.0.OR.IDIM( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX
            GO TO 220
         }
      } // 10
      // Values of ALPHA
      READ( NIN, FMT = * )NALF
      if ( NALF.LT.1.OR.NALF.GT.NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX
         GO TO 220
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF )
      // Values of BETA
      READ( NIN, FMT = * )NBET
      if ( NBET.LT.1.OR.NBET.GT.NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX
         GO TO 220
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET )

      // Report values of parameters.

      WRITE( NOUT, FMT = 9995 )
      WRITE( NOUT, FMT = 9994 )( IDIM( I ), I = 1, NIDIM )
      WRITE( NOUT, FMT = 9993 )( ALF( I ), I = 1, NALF )
      WRITE( NOUT, FMT = 9992 )( BET( I ), I = 1, NBET )
      if ( .NOT.TSTERR ) {
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9984 )
      }
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9999 )THRESH
      WRITE( NOUT, FMT = * )

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      for (I = 1; I <= NSUBS; I++) { // 20
         LTEST( I ) = .FALSE.
      } // 20
   30 READ( NIN, FMT = 9988, END = 60 )SNAMET, LTESTT
      for (I = 1; I <= NSUBS; I++) { // 40
         IF( SNAMET.EQ.SNAMES( I ) ) GO TO 50
      } // 40
      WRITE( NOUT, FMT = 9990 )SNAMET
      STOP
   50 LTEST( I ) = LTESTT
      GO TO 30

      } // 60
      CLOSE ( NIN )

      // Compute EPS (the machine precision).

      EPS = EPSILON(ZERO)
      WRITE( NOUT, FMT = 9998 )EPS

      // Check the reliability of DMMCH using exact data.

      N = MIN( 32, NMAX )
      for (J = 1; J <= N; J++) { // 100
         for (I = 1; I <= N; I++) { // 90
            AB( I, J ) = MAX( I - J + 1, 0 )
         } // 90
         AB( J, NMAX + 1 ) = J
         AB( 1, NMAX + J ) = J
         C( J, 1 ) = ZERO
      } // 100
      for (J = 1; J <= N; J++) { // 110
         CC( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
      } // 110
      // CC holds the exact result. On exit from DMMCH CT holds
      // the result computed by DMMCH.
      TRANSA = 'N'
      TRANSB = 'N'
      dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LDE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      TRANSB = 'T'
      dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LDE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      for (J = 1; J <= N; J++) { // 120
         AB( J, NMAX + 1 ) = N - J + 1
         AB( 1, NMAX + J ) = N - J + 1
      } // 120
      for (J = 1; J <= N; J++) { // 130
         CC( N - J + 1 ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
      } // 130
      TRANSA = 'T'
      TRANSB = 'N'
      dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LDE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      TRANSB = 'T'
      dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LDE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.ZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 200
         WRITE( NOUT, FMT = * )
         if ( .NOT.LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9987 )SNAMES( ISNUM )
         } else {
            SRNAMT = SNAMES( ISNUM )
            // Test error exits.
            if ( TSTERR ) {
               dchke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * )
            }
            // Test computations.
            INFOT = 0
            OK = .TRUE.
            FATAL = .FALSE.
            GO TO ( 140, 150, 160, 160, 170, 180 )ISNUM
            // Test DGEMM, 01.
  140       CALL DCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G )
            GO TO 190
            // Test DSYMM, 02.
  150       CALL DCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G )
            GO TO 190
            // Test DTRMM, 03, DTRSM, 04.
  160       CALL DCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C )
            GO TO 190
            // Test DSYRK, 05.
  170       CALL DCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G )
            GO TO 190
            // Test DSYR2K, 06.
  180       CALL DCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W )
            GO TO 190

  190       IF( FATAL.AND.SFATAL ) GO TO 210
         }
      } // 200
      WRITE( NOUT, FMT = 9986 )
      GO TO 230

      } // 210
      WRITE( NOUT, FMT = 9985 )
      GO TO 230

      } // 220
      WRITE( NOUT, FMT = 9991 )

      } // 230
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT )
      STOP

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 )
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, D9.1 )
 9997 FORMAT( ' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 )
 9995 FORMAT( ' TESTS OF THE double           LEVEL 3 BLAS', //' THE F',; 'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9994 FORMAT( '   FOR N              ', 9I6 )
 9993 FORMAT( '   FOR ALPHA          ', 7F6.1 )
 9992 FORMAT( '   FOR BETA           ', 7F6.1 )
 9991 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' )
 9990 FORMAT( ' SUBPROGRAM NAME ', A6, ' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' )
 9989 FORMAT( ' ERROR IN DMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' DMMCH WAS CALLED WITH TRANSA = ', A1, ' AND TRANSB = ', A1, /' AND RETURNED SAME = ', L1, ' AND ', 'ERR = ', F12.3, '.', /' THIS MAY BE DUE TO FAULTS IN THE ', 'ARITHMETIC OR THE COMPILER.', /' ******* TESTS ABANDONED ', '*******' )
 9988 FORMAT( A6, L2 )
 9987 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9986 FORMAT( /' END OF TESTS' )
 9985 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' )
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )

      // End of DBLAT3

      }
      SUBROUTINE DCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G )

*  Tests DGEMM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double             A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double             ALPHA, ALS, BETA, BLS, ERR, ERRMAX;
      int                I, IA, IB, ICA, ICB, IK, IM, IN, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MA, MB, MS, N, NA, NARGS, NB, NC, NS;
      bool               NULL, RESET, SAME, TRANA, TRANB;
      String             TRANAS, TRANBS, TRANSA, TRANSB;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LDE, LDERES;
      // EXTERNAL LDE, LDERES
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DMAKE, DMMCH
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'NTC'/
      // .. Executable Statements ..

      NARGS = 13
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO

      for (IM = 1; IM <= NIDIM; IM++) { // 110
         M = IDIM( IM )

         for (IN = 1; IN <= NIDIM; IN++) { // 100
            N = IDIM( IN )
            // Set LDC to 1 more than minimum value if room.
            LDC = M
            if (LDC.LT.NMAX) LDC = LDC + 1;
            // Skip tests if not enough room.
            if (LDC.GT.NMAX) GO TO 100;
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0

            for (IK = 1; IK <= NIDIM; IK++) { // 90
               K = IDIM( IK )

               for (ICA = 1; ICA <= 3; ICA++) { // 80
                  TRANSA = ICH( ICA: ICA )
                  TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'

                  if ( TRANA ) {
                     MA = K
                     NA = M
                  } else {
                     MA = M
                     NA = K
                  }
                  // Set LDA to 1 more than minimum value if room.
                  LDA = MA
                  if (LDA.LT.NMAX) LDA = LDA + 1;
                  // Skip tests if not enough room.
                  if (LDA.GT.NMAX) GO TO 80;
                  LAA = LDA*NA

                  // Generate the matrix A.

                  dmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (ICB = 1; ICB <= 3; ICB++) { // 70
                     TRANSB = ICH( ICB: ICB )
                     TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'

                     if ( TRANB ) {
                        MB = N
                        NB = K
                     } else {
                        MB = K
                        NB = N
                     }
                     // Set LDB to 1 more than minimum value if room.
                     LDB = MB
                     if (LDB.LT.NMAX) LDB = LDB + 1;
                     // Skip tests if not enough room.
                     if (LDB.GT.NMAX) GO TO 70;
                     LBB = LDB*NB

                     // Generate the matrix B.

                     dmake('GE', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO );

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA )

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB )

                           // Generate the matrix C.

                           dmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           TRANAS = TRANSA
                           TRANBS = TRANSB
                           MS = M
                           NS = N
                           KS = K
                           ALS = ALPHA
                           for (I = 1; I <= LAA; I++) { // 10
                              AS( I ) = AA( I )
                           } // 10
                           LDAS = LDA
                           for (I = 1; I <= LBB; I++) { // 20
                              BS( I ) = BB( I )
                           } // 20
                           LDBS = LDB
                           BLS = BETA
                           for (I = 1; I <= LCC; I++) { // 30
                              CS( I ) = CC( I )
                           } // 30
                           LDCS = LDC

                           // Call the subroutine.

                           if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC;
                           if (REWI) REWIND NTRA;
                           dgemm(TRANSA, TRANSB, M, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                           // Check if error-exit was taken incorrectly.

                           if ( .NOT.OK ) {
                              WRITE( NOUT, FMT = 9994 )
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = TRANSA.EQ.TRANAS
                           ISAME( 2 ) = TRANSB.EQ.TRANBS
                           ISAME( 3 ) = MS.EQ.M
                           ISAME( 4 ) = NS.EQ.N
                           ISAME( 5 ) = KS.EQ.K
                           ISAME( 6 ) = ALS.EQ.ALPHA
                           ISAME( 7 ) = LDE( AS, AA, LAA )
                           ISAME( 8 ) = LDAS.EQ.LDA
                           ISAME( 9 ) = LDE( BS, BB, LBB )
                           ISAME( 10 ) = LDBS.EQ.LDB
                           ISAME( 11 ) = BLS.EQ.BETA
                           if ( NULL ) {
                              ISAME( 12 ) = LDE( CS, CC, LCC )
                           } else {
                              ISAME( 12 ) = LDERES( 'GE', ' ', M, N, CS, CC, LDC )
                           }
                           ISAME( 13 ) = LDCS.EQ.LDC

                           // If data was incorrectly changed, report
                           // and return.

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

                              dmmch(TRANSA, TRANSB, M, N, K, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                              ERRMAX = MAX( ERRMAX, ERR )
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC

      } // 130
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',''', A1, ''',', 3( I3, ',' ), F4.1, ', A,', I3, ', B,', I3, ',', F4.1, ', ', 'C,', I3, ').' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of DCHK1

      }
      SUBROUTINE DCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G )

*  Tests DSYMM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double             A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double             ALPHA, ALS, BETA, BLS, ERR, ERRMAX;
      int                I, IA, IB, ICS, ICU, IM, IN, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String             SIDE, SIDES, UPLO, UPLOS;
      String             ICHS, ICHU;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LDE, LDERES;
      // EXTERNAL LDE, LDERES
      // .. External Subroutines ..
      // EXTERNAL DMAKE, DMMCH, DSYMM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHS/'LR'/, ICHU/'UL'/
      // .. Executable Statements ..

      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO

      for (IM = 1; IM <= NIDIM; IM++) { // 100
         M = IDIM( IM )

         for (IN = 1; IN <= NIDIM; IN++) { // 90
            N = IDIM( IN )
            // Set LDC to 1 more than minimum value if room.
            LDC = M
            if (LDC.LT.NMAX) LDC = LDC + 1;
            // Skip tests if not enough room.
            if (LDC.GT.NMAX) GO TO 90;
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0

            // Set LDB to 1 more than minimum value if room.
            LDB = M
            if (LDB.LT.NMAX) LDB = LDB + 1;
            // Skip tests if not enough room.
            if (LDB.GT.NMAX) GO TO 90;
            LBB = LDB*N

            // Generate the matrix B.

            dmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

            for (ICS = 1; ICS <= 2; ICS++) { // 80
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'

               if ( LEFT ) {
                  NA = M
               } else {
                  NA = N
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA
               if (LDA.LT.NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA.GT.NMAX) GO TO 80;
               LAA = LDA*NA

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU )

                  // Generate the symmetric matrix A.

                  dmake('SY', UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA )

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB )

                        // Generate the matrix C.

                        dmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1

                        // Save every datum before calling the
                        // subroutine.

                        SIDES = SIDE
                        UPLOS = UPLO
                        MS = M
                        NS = N
                        ALS = ALPHA
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I )
                        } // 10
                        LDAS = LDA
                        for (I = 1; I <= LBB; I++) { // 20
                           BS( I ) = BB( I )
                        } // 20
                        LDBS = LDB
                        BLS = BETA
                        for (I = 1; I <= LCC; I++) { // 30
                           CS( I ) = CC( I )
                        } // 30
                        LDCS = LDC

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        dsymm(SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9994 )
                           FATAL = .TRUE.
                           GO TO 110
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = SIDES.EQ.SIDE
                        ISAME( 2 ) = UPLOS.EQ.UPLO
                        ISAME( 3 ) = MS.EQ.M
                        ISAME( 4 ) = NS.EQ.N
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LDE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LDE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        ISAME( 10 ) = BLS.EQ.BETA
                        if ( NULL ) {
                           ISAME( 11 ) = LDE( CS, CC, LCC )
                        } else {
                           ISAME( 11 ) = LDERES( 'GE', ' ', M, N, CS, CC, LDC )
                        }
                        ISAME( 12 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                        } // 40
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 110
                        }

                        if ( .NOT.NULL ) {

                           // Check the result.

                           if ( LEFT ) {
                              dmmch('N', 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                           } else {
                              dmmch('N', 'N', M, N, N, ALPHA, B, NMAX, A, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 120

      } // 110
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC

      } // 120
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ), F4.1, ', A,', I3, ', B,', I3, ',', F4.1, ', C,', I3, ')   ', ' .' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of DCHK2

      }
      SUBROUTINE DCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, A, AA, AS, B, BB, BS, CT, G, C )

*  Tests DTRMM and DTRSM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                NALF, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double             A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double             ALPHA, ALS, ERR, ERRMAX;
      int                I, IA, ICD, ICS, ICT, ICU, IM, IN, J, LAA, LBB, LDA, LDAS, LDB, LDBS, M, MS, N, NA, NARGS, NC, NS;
      bool               LEFT, NULL, RESET, SAME;
      String             DIAG, DIAGS, SIDE, SIDES, TRANAS, TRANSA, UPLO, UPLOS;
      String             ICHD, ICHS, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LDE, LDERES;
      // EXTERNAL LDE, LDERES
      // .. External Subroutines ..
      // EXTERNAL DMAKE, DMMCH, DTRMM, DTRSM
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/, ICHS/'LR'/
      // .. Executable Statements ..

      NARGS = 11
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
      // Set up zero matrix for DMMCH.
      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            C( I, J ) = ZERO
         } // 10
      } // 20

      for (IM = 1; IM <= NIDIM; IM++) { // 140
         M = IDIM( IM )

         for (IN = 1; IN <= NIDIM; IN++) { // 130
            N = IDIM( IN )
            // Set LDB to 1 more than minimum value if room.
            LDB = M
            if (LDB.LT.NMAX) LDB = LDB + 1;
            // Skip tests if not enough room.
            if (LDB.GT.NMAX) GO TO 130;
            LBB = LDB*N
            NULL = M.LE.0.OR.N.LE.0

            for (ICS = 1; ICS <= 2; ICS++) { // 120
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'
               if ( LEFT ) {
                  NA = M
               } else {
                  NA = N
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA
               if (LDA.LT.NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA.GT.NMAX) GO TO 130;
               LAA = LDA*NA

               for (ICU = 1; ICU <= 2; ICU++) { // 110
                  UPLO = ICHU( ICU: ICU )

                  for (ICT = 1; ICT <= 3; ICT++) { // 100
                     TRANSA = ICHT( ICT: ICT )

                     for (ICD = 1; ICD <= 2; ICD++) { // 90
                        DIAG = ICHD( ICD: ICD )

                        for (IA = 1; IA <= NALF; IA++) { // 80
                           ALPHA = ALF( IA )

                           // Generate the matrix A.

                           dmake('TR', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO );

                           // Generate the matrix B.

                           dmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO );

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           SIDES = SIDE
                           UPLOS = UPLO
                           TRANAS = TRANSA
                           DIAGS = DIAG
                           MS = M
                           NS = N
                           ALS = ALPHA
                           for (I = 1; I <= LAA; I++) { // 30
                              AS( I ) = AA( I )
                           } // 30
                           LDAS = LDA
                           for (I = 1; I <= LBB; I++) { // 40
                              BS( I ) = BB( I )
                           } // 40
                           LDBS = LDB

                           // Call the subroutine.

                           if ( SNAME( 4: 5 ).EQ.'MM' ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;
                              if (REWI) REWIND NTRA;
                              dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           } else if ( SNAME( 4: 5 ).EQ.'SM' ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB;
                              if (REWI) REWIND NTRA;
                              dtrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( .NOT.OK ) {
                              WRITE( NOUT, FMT = 9994 )
                              FATAL = .TRUE.
                              GO TO 150
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = SIDES.EQ.SIDE
                           ISAME( 2 ) = UPLOS.EQ.UPLO
                           ISAME( 3 ) = TRANAS.EQ.TRANSA
                           ISAME( 4 ) = DIAGS.EQ.DIAG
                           ISAME( 5 ) = MS.EQ.M
                           ISAME( 6 ) = NS.EQ.N
                           ISAME( 7 ) = ALS.EQ.ALPHA
                           ISAME( 8 ) = LDE( AS, AA, LAA )
                           ISAME( 9 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 10 ) = LDE( BS, BB, LBB )
                           } else {
                              ISAME( 10 ) = LDERES( 'GE', ' ', M, N, BS, BB, LDB )
                           }
                           ISAME( 11 ) = LDBS.EQ.LDB

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = .TRUE.
                           for (I = 1; I <= NARGS; I++) { // 50
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                           } // 50
                           if ( .NOT.SAME ) {
                              FATAL = .TRUE.
                              GO TO 150
                           }

                           if ( .NOT.NULL ) {
                              if ( SNAME( 4: 5 ).EQ.'MM' ) {

                                 // Check the result.

                                 if ( LEFT ) {
                                    dmmch(TRANSA, 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .TRUE. );
                                 } else {
                                    dmmch('N', TRANSA, M, N, N, ALPHA, B, NMAX, A, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .TRUE. );
                                 }
                              } else if ( SNAME( 4: 5 ).EQ.'SM' ) {

                                 // Compute approximation to original
                                 // matrix.

                                 for (J = 1; J <= N; J++) { // 70
                                    for (I = 1; I <= M; I++) { // 60
                                       C( I, J ) = BB( I + ( J - 1 )* LDB )                                        BB( I + ( J - 1 )*LDB ) = ALPHA* B( I, J )
                                    } // 60
                                 } // 70

                                 if ( LEFT ) {
                                    dmmch(TRANSA, 'N', M, N, M, ONE, A, NMAX, C, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .FALSE. );
                                 } else {
                                    dmmch('N', TRANSA, M, N, N, ONE, C, NMAX, A, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .FALSE. );
                                 }
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 160

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB

      } // 160
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(', 4( '''', A1, ''',' ), 2( I3, ',' ), F4.1, ', A,', I3, ', B,', I3, ')        .' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of DCHK3

      }
      SUBROUTINE DCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G )

*  Tests DSYRK.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double             A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double             ALPHA, ALS, BETA, BETS, ERR, ERRMAX;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, KS, LAA, LCC, LDA, LDAS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, UPLO, UPLOS;
      String             ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LDE, LDERES;
      // EXTERNAL LDE, LDERES
      // .. External Subroutines ..
      // EXTERNAL DMAKE, DMMCH, DSYRK
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NTC'/, ICHU/'UL'/
      // .. Executable Statements ..

      NARGS = 10
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN )
         // Set LDC to 1 more than minimum value if room.
         LDC = N
         if (LDC.LT.NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC.GT.NMAX) GO TO 100;
         LCC = LDC*N
         NULL = N.LE.0

         for (IK = 1; IK <= NIDIM; IK++) { // 90
            K = IDIM( IK )

            for (ICT = 1; ICT <= 3; ICT++) { // 80
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'
               if ( TRAN ) {
                  MA = K
                  NA = N
               } else {
                  MA = N
                  NA = K
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA
               if (LDA.LT.NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA.GT.NMAX) GO TO 80;
               LAA = LDA*NA

               // Generate the matrix A.

               dmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO );

               for (ICU = 1; ICU <= 2; ICU++) { // 70
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'

                  for (IA = 1; IA <= NALF; IA++) { // 60
                     ALPHA = ALF( IA )

                     for (IB = 1; IB <= NBET; IB++) { // 50
                        BETA = BET( IB )

                        // Generate the matrix C.

                        dmake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        ALS = ALPHA
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I )
                        } // 10
                        LDAS = LDA
                        BETS = BETA
                        for (I = 1; I <= LCC; I++) { // 20
                           CS( I ) = CC( I )
                        } // 20
                        LDCS = LDC

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        dsyrk(UPLO, TRANS, N, K, ALPHA, AA, LDA, BETA, CC, LDC );

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9993 )
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS.EQ.UPLO
                        ISAME( 2 ) = TRANSS.EQ.TRANS
                        ISAME( 3 ) = NS.EQ.N
                        ISAME( 4 ) = KS.EQ.K
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LDE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = BETS.EQ.BETA
                        if ( NULL ) {
                           ISAME( 9 ) = LDE( CS, CC, LCC )
                        } else {
                           ISAME( 9 ) = LDERES( 'SY', UPLO, N, N, CS, CC, LDC )
                        }
                        ISAME( 10 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

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

                           JC = 1
                           for (J = 1; J <= N; J++) { // 40
                              if ( UPPER ) {
                                 JJ = 1
                                 LJ = J
                              } else {
                                 JJ = J
                                 LJ = N - J + 1
                              }
                              if ( TRAN ) {
                                 dmmch('T', 'N', LJ, 1, K, ALPHA, A( 1, JJ ), NMAX, A( 1, J ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                              } else {
                                 dmmch('N', 'T', LJ, 1, K, ALPHA, A( JJ, 1 ), NMAX, A( J, 1 ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC
                              } else {
                                 JC = JC + LDC + 1
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 110
      if (N.GT.1) WRITE( NOUT, FMT = 9995 )J;

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC

      } // 130
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ), F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ')           .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of DCHK4

      }
      SUBROUTINE DCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W )

*  Tests DSYR2K.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      double             AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), G( NMAX ), W( 2*NMAX );
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      double             ALPHA, ALS, BETA, BETS, ERR, ERRMAX;
      int                I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB, K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NS;
      bool               NULL, RESET, SAME, TRAN, UPPER;
      String             TRANS, TRANSS, UPLO, UPLOS;
      String             ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LDE, LDERES;
      // EXTERNAL LDE, LDERES
      // .. External Subroutines ..
      // EXTERNAL DMAKE, DMMCH, DSYR2K
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NTC'/, ICHU/'UL'/
      // .. Executable Statements ..

      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 130
         N = IDIM( IN )
         // Set LDC to 1 more than minimum value if room.
         LDC = N
         if (LDC.LT.NMAX) LDC = LDC + 1;
         // Skip tests if not enough room.
         if (LDC.GT.NMAX) GO TO 130;
         LCC = LDC*N
         NULL = N.LE.0

         for (IK = 1; IK <= NIDIM; IK++) { // 120
            K = IDIM( IK )

            for (ICT = 1; ICT <= 3; ICT++) { // 110
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'
               if ( TRAN ) {
                  MA = K
                  NA = N
               } else {
                  MA = N
                  NA = K
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA
               if (LDA.LT.NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA.GT.NMAX) GO TO 110;
               LAA = LDA*NA

               // Generate the matrix A.

               if ( TRAN ) {
                  dmake('GE', ' ', ' ', MA, NA, AB, 2*NMAX, AA, LDA, RESET, ZERO );
               } else {
                  dmake('GE', ' ', ' ', MA, NA, AB, NMAX, AA, LDA, RESET, ZERO );
               }

               // Generate the matrix B.

               LDB = LDA
               LBB = LAA
               if ( TRAN ) {
                  dmake('GE', ' ', ' ', MA, NA, AB( K + 1 ), 2*NMAX, BB, LDB, RESET, ZERO );
               } else {
                  dmake('GE', ' ', ' ', MA, NA, AB( K*NMAX + 1 ), NMAX, BB, LDB, RESET, ZERO );
               }

               for (ICU = 1; ICU <= 2; ICU++) { // 100
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'

                  for (IA = 1; IA <= NALF; IA++) { // 90
                     ALPHA = ALF( IA )

                     for (IB = 1; IB <= NBET; IB++) { // 80
                        BETA = BET( IB )

                        // Generate the matrix C.

                        dmake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO );

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        ALS = ALPHA
                        for (I = 1; I <= LAA; I++) { // 10
                           AS( I ) = AA( I )
                        } // 10
                        LDAS = LDA
                        for (I = 1; I <= LBB; I++) { // 20
                           BS( I ) = BB( I )
                        } // 20
                        LDBS = LDB
                        BETS = BETA
                        for (I = 1; I <= LCC; I++) { // 30
                           CS( I ) = CC( I )
                        } // 30
                        LDCS = LDC

                        // Call the subroutine.

                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC;
                        if (REWI) REWIND NTRA;
                        dsyr2k(UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC );

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9993 )
                           FATAL = .TRUE.
                           GO TO 150
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS.EQ.UPLO
                        ISAME( 2 ) = TRANSS.EQ.TRANS
                        ISAME( 3 ) = NS.EQ.N
                        ISAME( 4 ) = KS.EQ.K
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LDE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LDE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        ISAME( 10 ) = BETS.EQ.BETA
                        if ( NULL ) {
                           ISAME( 11 ) = LDE( CS, CC, LCC )
                        } else {
                           ISAME( 11 ) = LDERES( 'SY', UPLO, N, N, CS, CC, LDC )
                        }
                        ISAME( 12 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                        } // 40
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 150
                        }

                        if ( .NOT.NULL ) {

                           // Check the result column by column.

                           JJAB = 1
                           JC = 1
                           for (J = 1; J <= N; J++) { // 70
                              if ( UPPER ) {
                                 JJ = 1
                                 LJ = J
                              } else {
                                 JJ = J
                                 LJ = N - J + 1
                              }
                              if ( TRAN ) {
                                 for (I = 1; I <= K; I++) { // 50
                                    W( I ) = AB( ( J - 1 )*2*NMAX + K + I )                                     W( K + I ) = AB( ( J - 1 )*2*NMAX + I )
                                 } // 50
                                 dmmch('T', 'N', LJ, 1, 2*K, ALPHA, AB( JJAB ), 2*NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                              } else {
                                 for (I = 1; I <= K; I++) { // 60
                                    W( I ) = AB( ( K + I - 1 )*NMAX + J )                                     W( K + I ) = AB( ( I - 1 )*NMAX + J )
                                 } // 60
                                 dmmch('N', 'N', LJ, 1, 2*K, ALPHA, AB( JJ ), NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. );
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC
                              } else {
                                 JC = JC + LDC + 1
                                 if (TRAN) JJAB = JJAB + 2*NMAX;
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
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

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 160

      } // 140
      if (N.GT.1) WRITE( NOUT, FMT = 9995 )J;

      } // 150
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC

      } // 160
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ), F4.1, ', A,', I3, ', B,', I3, ',', F4.1, ', C,', I3, ')   ', ' .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of DCHK5

      }
      SUBROUTINE DCHKE( ISNUM, SRNAMT, NOUT )

*  Tests the error exits from the Level 3 Blas.
*  Requires a special version of the error-handling routine XERBLA.
*  A, B and C should not need to be defined.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

*  3-19-92:  Initialize ALPHA and BETA  (eca)
*  3-19-92:  Fix argument 12 in calls to SSYMM with INFOT = 9  (eca)

      // .. Scalar Arguments ..
      int                ISNUM, NOUT;
      String             SRNAMT;
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Parameters ..
      double             ONE, TWO;
      const              ONE = 1.0D0, TWO = 2.0D0 ;
      // .. Local Scalars ..
      double             ALPHA, BETA;
      // .. Local Arrays ..
      double             A( 2, 1 ), B( 2, 1 ), C( 2, 1 );
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGEMM, DSYMM, DSYR2K, DSYRK, DTRMM, DTRSM
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      // OK is set to .FALSE. by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      OK = .TRUE.
      // LERR is set to .TRUE. by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      LERR = .FALSE.

      // Initialize ALPHA and BETA.

      ALPHA = ONE
      BETA = TWO

      GO TO ( 10, 20, 30, 40, 50, 60 )ISNUM
   10 INFOT = 1
      dgemm('/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 1
      dgemm('/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgemm('N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgemm('T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgemm('N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgemm('N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgemm('T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgemm('T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgemm('N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgemm('N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgemm('T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgemm('T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemm('N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemm('N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemm('T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemm('T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemm('N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemm('N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemm('T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemm('T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgemm('N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgemm('T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgemm('N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgemm('T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      dgemm('N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      dgemm('N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      dgemm('T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      dgemm('T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70
   20 INFOT = 1
      dsymm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dsymm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsymm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsymm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsymm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsymm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsymm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsymm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsymm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsymm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsymm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsymm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsymm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsymm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70
   30 INFOT = 1
      dtrmm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrmm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dtrmm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrmm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrmm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrmm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrmm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrmm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70
   40 INFOT = 1
      dtrsm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dtrsm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dtrsm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dtrsm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      dtrsm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      dtrsm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dtrsm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      dtrsm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70
   50 INFOT = 1
      dsyrk('/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dsyrk('U', '/', 0, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyrk('U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyrk('U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyrk('L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyrk('L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyrk('U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyrk('U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyrk('L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyrk('L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyrk('U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyrk('U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyrk('L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyrk('L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dsyrk('U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dsyrk('U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dsyrk('L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      dsyrk('L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 70
   60 INFOT = 1
      dsyr2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      dsyr2k('U', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyr2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyr2k('U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyr2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      dsyr2k('L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyr2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyr2k('U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyr2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      dsyr2k('L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyr2k('U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      dsyr2k('L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsyr2k('U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      dsyr2k('L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsyr2k('U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 12
      dsyr2k('L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );

   70 IF( OK )THEN
         WRITE( NOUT, FMT = 9999 )SRNAMT
      } else {
         WRITE( NOUT, FMT = 9998 )SRNAMT
      }
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE TESTS OF ERROR-EXITS' )
 9998 FORMAT( ' ******* ', A6, ' FAILED THE TESTS OF ERROR-EXITS *****', '**' )

      // End of DCHKE

      }
      SUBROUTINE DMAKE( TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, RESET, TRANSL )

*  Generates values for an M by N matrix A.
*  Stores the values in the array AA in the data structure required
*  by the routine, with unwanted elements set to rogue value.

*  TYPE is 'GE', 'SY' or 'TR'.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      double             ROGUE;
      const              ROGUE = -1.0D10 ;
      // .. Scalar Arguments ..
      double             TRANSL;
      int                LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      double             A( NMAX, * ), AA( * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      double             DBEG;
      // EXTERNAL DBEG
      // .. Executable Statements ..
      GEN = TYPE.EQ.'GE'
      SYM = TYPE.EQ.'SY'
      TRI = TYPE.EQ.'TR'
      UPPER = ( SYM.OR.TRI ).AND.UPLO.EQ.'U'
      LOWER = ( SYM.OR.TRI ).AND.UPLO.EQ.'L'
      UNIT = TRI.AND.DIAG.EQ.'U'

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) ) {
               A( I, J ) = DBEG( RESET ) + TRANSL
               if ( I.NE.J ) {
                  // Set some elements to zero
                  if (N.GT.3.AND.J.EQ.N/2) A( I, J ) = ZERO;
                  if ( SYM ) {
                     A( J, I ) = A( I, J )
                  } else if ( TRI ) {
                     A( J, I ) = ZERO
                  }
               }
            }
         } // 10
         if (TRI) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE;
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE.EQ.'GE' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 30
            for (I = M + 1; I <= LDA; I++) { // 40
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 40
         } // 50
      } else if ( TYPE.EQ.'SY'.OR.TYPE.EQ.'TR' ) {
         for (J = 1; J <= N; J++) { // 90
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
            for (I = 1; I <= IBEG - 1; I++) { // 60
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 60
            for (I = IBEG; I <= IEND; I++) { // 70
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 70
            for (I = IEND + 1; I <= LDA; I++) { // 80
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 80
         } // 90
      }
      RETURN

      // End of DMAKE

      }
      SUBROUTINE DMMCH( TRANSA, TRANSB, M, N, KK, ALPHA, A, LDA, B, LDB, BETA, C, LDC, CT, G, CC, LDCC, EPS, ERR, FATAL, NOUT, MV )

*  Checks the results of the computational tests.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // .. Scalar Arguments ..
      double             ALPHA, BETA, EPS, ERR;
      int                KK, LDA, LDB, LDC, LDCC, M, N, NOUT;
      bool               FATAL, MV;
      String             TRANSA, TRANSB;
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), C( LDC, * ), CC( LDCC, * ), CT( * ), G( * );
      // .. Local Scalars ..
      double             ERRI;
      int                I, J, K;
      bool               TRANA, TRANB;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // .. Executable Statements ..
      TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'
      TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'

      // Compute expected result, one column at a time, in CT using data
      // in A, B and C.
      // Compute gauges in G.

      for (J = 1; J <= N; J++) { // 120

         for (I = 1; I <= M; I++) { // 10
            CT( I ) = ZERO
            G( I ) = ZERO
         } // 10
         if ( .NOT.TRANA.AND..NOT.TRANB ) {
            for (K = 1; K <= KK; K++) { // 30
               for (I = 1; I <= M; I++) { // 20
                  CT( I ) = CT( I ) + A( I, K )*B( K, J )
                  G( I ) = G( I ) + ABS( A( I, K ) )*ABS( B( K, J ) )
               } // 20
            } // 30
         } else if ( TRANA.AND..NOT.TRANB ) {
            for (K = 1; K <= KK; K++) { // 50
               for (I = 1; I <= M; I++) { // 40
                  CT( I ) = CT( I ) + A( K, I )*B( K, J )
                  G( I ) = G( I ) + ABS( A( K, I ) )*ABS( B( K, J ) )
               } // 40
            } // 50
         } else if ( .NOT.TRANA.AND.TRANB ) {
            for (K = 1; K <= KK; K++) { // 70
               for (I = 1; I <= M; I++) { // 60
                  CT( I ) = CT( I ) + A( I, K )*B( J, K )
                  G( I ) = G( I ) + ABS( A( I, K ) )*ABS( B( J, K ) )
               } // 60
            } // 70
         } else if ( TRANA.AND.TRANB ) {
            for (K = 1; K <= KK; K++) { // 90
               for (I = 1; I <= M; I++) { // 80
                  CT( I ) = CT( I ) + A( K, I )*B( J, K )
                  G( I ) = G( I ) + ABS( A( K, I ) )*ABS( B( J, K ) )
               } // 80
            } // 90
         }
         for (I = 1; I <= M; I++) { // 100
            CT( I ) = ALPHA*CT( I ) + BETA*C( I, J )
            G( I ) = ABS( ALPHA )*G( I ) + ABS( BETA )*ABS( C( I, J ) )
         } // 100

         // Compute the error ratio for this result.

         ERR = ZERO
         for (I = 1; I <= M; I++) { // 110
            ERRI = ABS( CT( I ) - CC( I, J ) )/EPS
            IF( G( I ).NE.ZERO ) ERRI = ERRI/G( I )
            ERR = MAX( ERR, ERRI )
            IF( ERR*SQRT( EPS ).GE.ONE ) GO TO 130
         } // 110

      } // 120

      // If the loop completes, all results are at least half accurate.
      GO TO 150

      // Report fatal error.

  130 FATAL = .TRUE.
      WRITE( NOUT, FMT = 9999 )
      for (I = 1; I <= M; I++) { // 140
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, CT( I ), CC( I, J )
         } else {
            WRITE( NOUT, FMT = 9998 )I, CC( I, J ), CT( I )
         }
      } // 140
      if (N.GT.1) WRITE( NOUT, FMT = 9997 )J;

      } // 150
      RETURN

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'           EXPECTED RESULT   COMPU', 'TED RESULT' )
 9998 FORMAT( 1X, I7, 2G18.6 )
 9997 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )

      // End of DMMCH

      }
      bool    FUNCTION LDE( RI, RJ, LR );

*  Tests if two arrays are identical.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      double             RI( * ), RJ( * );
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         IF( RI( I ).NE.RJ( I ) ) GO TO 20
      } // 10
      LDE = .TRUE.
      GO TO 30
      } // 20
      LDE = .FALSE.
   30 RETURN

      // End of LDE

      }
      bool    FUNCTION LDERES( TYPE, UPLO, M, N, AA, AS, LDA );

*  Tests if selected elements in two arrays are equal.

*  TYPE is 'GE' or 'SY'.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      double             AA( LDA, * ), AS( LDA, * );
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO.EQ.'U'
      if ( TYPE.EQ.'GE' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 10
         } // 20
      } else if ( TYPE.EQ.'SY' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            for (I = 1; I <= IBEG - 1; I++) { // 30
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 30
            for (I = IEND + 1; I <= LDA; I++) { // 40
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
            } // 40
         } // 50
      }

      LDERES = .TRUE.
      GO TO 80
      } // 70
      LDERES = .FALSE.
   80 RETURN

      // End of LDERES

      }
      double           FUNCTION DBEG( RESET );

*  Generates random numbers uniformly distributed between -0.5 and 0.5.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, MI;
      // .. Save statement ..
      SAVE               I, IC, MI
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891
         I = 7
         IC = 0
         RESET = .FALSE.
      }

      // The sequence of values of I is bounded between 1 and 999.
      // If initial I = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I = 4 or 8, the period will be 25.
      // If initial I = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I in 6.

      IC = IC + 1
   10 I = I*MI
      I = I - 1000*( I/1000 )
      if ( IC.GE.5 ) {
         IC = 0
         GO TO 10
      }
      DBEG = ( I - 500 )/1001.0D0
      RETURN

      // End of DBEG

      }
      double           FUNCTION DDIFF( X, Y );

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      double             X, Y;
      // .. Executable Statements ..
      DDIFF = X - Y
      RETURN

      // End of DDIFF

      }
      SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK )

*  Tests whether XERBLA has detected an error when it should.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Executable Statements ..
      if ( .NOT.LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT, SRNAMT
         OK = .FALSE.
      }
      LERR = .FALSE.
      RETURN

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ', I2, ' NOT D', 'ETECTED BY ', A6, ' *****' )

      // End of CHKXER

      }
      SUBROUTINE XERBLA( SRNAME, INFO )

*  This is a special version of XERBLA to be used only as part of
*  the test program for testing error exits from the Level 3 BLAS
*  routines.

*  XERBLA  is an error handler for the Level 3 BLAS routines.

*  It is called by the Level 3 BLAS routines if an input parameter is
*  invalid.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      int                INFO;
      String             SRNAME;
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUT, OK, LERR
      // COMMON /SRNAMC/SRNAMT
      // .. Executable Statements ..
      LERR = .TRUE.
      if ( INFO.NE.INFOT ) {
         if ( INFOT.NE.0 ) {
            WRITE( NOUT, FMT = 9999 )INFO, INFOT
         } else {
            WRITE( NOUT, FMT = 9997 )INFO
         }
         OK = .FALSE.
      }
      if ( SRNAME.NE.SRNAMT ) {
         WRITE( NOUT, FMT = 9998 )SRNAME, SRNAMT
         OK = .FALSE.
      }
      RETURN

 9999 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' INSTEAD', ' OF ', I2, ' *******' )
 9998 FORMAT( ' ******* XERBLA WAS CALLED WITH SRNAME = ', A6, ' INSTE', 'AD OF ', A6, ' *******' )
 9997 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' *******' )

      // End of XERBLA

      }
