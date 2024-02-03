      PROGRAM CBLAT3

*  Test program for the COMPLEX          Level 3 Blas.

*  The program must be driven by a short data file. The first 13 records
*  of the file are read using list-directed input, the last 9 records
*  are read using the format ( A12, L2 ). An annotated example of a data
*  file can be obtained by deleting the first 3 characters from the
*  following 22 lines:
*  'CBLAT3.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
*  -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
*  F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
*  F        LOGICAL FLAG, T TO STOP ON FAILURES.
*  T        LOGICAL FLAG, T TO TEST ERROR EXITS.
*  2        0 TO TEST COLUMN-MAJOR, 1 TO TEST ROW-MAJOR, 2 TO TEST BOTH
*  16.0     THRESHOLD VALUE OF TEST RATIO
*  6                 NUMBER OF VALUES OF N
*  0 1 2 3 5 9       VALUES OF N
*  3                 NUMBER OF VALUES OF ALPHA
*  (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
*  3                 NUMBER OF VALUES OF BETA
*  (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
*  cblas_cgemm  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_chemm  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_csymm  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ctrmm  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_ctrsm  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_cherk  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_csyrk  T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_cher2k T PUT F FOR NO TEST. SAME COLUMNS.
*  cblas_csyr2k T PUT F FOR NO TEST. SAME COLUMNS.

*  See:

      // Dongarra J. J., Du Croz J. J., Duff I. S. and Hammarling S.
      // A Set of Level 3 Basic Linear Algebra Subprograms.

      // Technical Memorandum No.88 (Revision 1), Mathematics and
      // Computer Science Division, Argonne National Laboratory, 9700
      // South Cass Avenue, Argonne, Illinois 60439, US.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      int                NIN, NOUT;
      const              NIN = 5, NOUT = 6 ;
      int                NSUBS;
      const              NSUBS = 9 ;
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO, RHALF, RONE
      const              RZERO = 0.0, RHALF = 0.5, RONE = 1.0 ;
      int                NMAX;
      const              NMAX = 65 ;
      int                NIDMAX, NALMAX, NBEMAX;
      const              NIDMAX = 9, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      REAL               EPS, ERR, THRESH
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NTRA, LAYOUT;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR, CORDER, RORDER;
      String             TRANSA, TRANSB;
      String             SNAMET;
      String             SNAPS;
      // .. Local Arrays ..
      COMPLEX            AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBEMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      REAL               SDIFF
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      COMMON             /SRNAMC/SRNAMT
      // .. Data statements ..
      DATA               SNAMES/'cblas_cgemm ', 'cblas_chemm ', 'cblas_csymm ', 'cblas_ctrmm ', 'cblas_ctrsm ', 'cblas_cherk ', 'cblas_csyrk ', 'cblas_cher2k', 'cblas_csyr2k'/
      // .. Executable Statements ..

      NOUTC = NOUT

      // Read name and unit number for snapshot output file and open file.

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
         GO TO 220
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM )
      DO 10 I = 1, NIDIM
         if ( IDIM( I ).LT.0.OR.IDIM( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX
            GO TO 220
         }
   10 CONTINUE
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

      DO 20 I = 1, NSUBS
         LTEST( I ) = .FALSE.
   20 CONTINUE
   30 READ( NIN, FMT = 9988, END = 60 )SNAMET, LTESTT
      DO 40 I = 1, NSUBS
         IF( SNAMET.EQ.SNAMES( I ) ) GO TO 50
   40 CONTINUE
      WRITE( NOUT, FMT = 9990 )SNAMET
      STOP
   50 LTEST( I ) = LTESTT
      GO TO 30

   60 CONTINUE
      CLOSE ( NIN )

      // Compute EPS (the machine precision).

      EPS = RONE
   70 CONTINUE
      IF( SDIFF( RONE + EPS, RONE ).EQ.RZERO ) GO TO 80
      EPS = RHALF*EPS
      GO TO 70
   80 CONTINUE
      EPS = EPS + EPS
      WRITE( NOUT, FMT = 9998 )EPS

      // Check the reliability of CMMCH using exact data.

      N = MIN( 32, NMAX )
      DO 100 J = 1, N
         DO 90 I = 1, N
            AB( I, J ) = MAX( I - J + 1, 0 )
   90    CONTINUE
         AB( J, NMAX + 1 ) = J
         AB( 1, NMAX + J ) = J
         C( J, 1 ) = ZERO
  100 CONTINUE
      DO 110 J = 1, N
         CC( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
  110 CONTINUE
      // CC holds the exact result. On exit from CMMCH CT holds
     t // he result computed by CMMCH.
      TRANSA = 'N'
      TRANSB = 'N'
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. )
      SAME = LCE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      TRANSB = 'C'
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. )
      SAME = LCE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      DO 120 J = 1, N
         AB( J, NMAX + 1 ) = N - J + 1
         AB( 1, NMAX + J ) = N - J + 1
  120 CONTINUE
      DO 130 J = 1, N
         CC( N - J + 1 ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
  130 CONTINUE
      TRANSA = 'C'
      TRANSB = 'N'
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. )
      SAME = LCE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }
      TRANSB = 'C'
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC, NMAX, EPS, ERR, FATAL, NOUT, .TRUE. )
      SAME = LCE( CC, CT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
         STOP
      }

      // Test each subroutine in turn.

      DO 200 ISNUM = 1, NSUBS
         WRITE( NOUT, FMT = * )
         if ( .NOT.LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9987 )SNAMES( ISNUM )
         } else {
            SRNAMT = SNAMES( ISNUM )
            // Test error exits.
            if ( TSTERR ) {
               CALL CC3CHKE( SNAMES( ISNUM ) )
               WRITE( NOUT, FMT = * )
            }
            // Test computations.
            INFOT = 0
            OK = .TRUE.
            FATAL = .FALSE.
            GO TO ( 140, 150, 150, 160, 160, 170, 170, 180, 180 )ISNUM
            // Test CGEMM, 01.
  140       IF (CORDER) THEN
            CALL CCHK1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 )
            }
            if (RORDER) {
            CALL CCHK1(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 )
            }
            GO TO 190
            // Test CHEMM, 02, CSYMM, 03.
  150       IF (CORDER) THEN
            CALL CCHK2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 )
            }
            if (RORDER) {
            CALL CCHK2(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 )
            }
            GO TO 190
            // Test CTRMM, 04, CTRSM, 05.
  160       IF (CORDER) THEN
            CALL CCHK3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C, 0 )
            }
            if (RORDER) {
            CALL CCHK3(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C, 1 )
            }
            GO TO 190
            // Test CHERK, 06, CSYRK, 07.
  170       IF (CORDER) THEN
            CALL CCHK4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 0 )
            }
            if (RORDER) {
            CALL CCHK4(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C, CC, CS, CT, G, 1 )
            }
            GO TO 190
            // Test CHER2K, 08, CSYR2K, 09.
  180       IF (CORDER) THEN
            CALL CCHK5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, 0 )
            }
            if (RORDER) {
            CALL CCHK5(SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, 1 )
            }
            GO TO 190

  190       IF( FATAL.AND.SFATAL )
     $         GO TO 210
         }
  200 CONTINUE
      WRITE( NOUT, FMT = 9986 )
      GO TO 230

  210 CONTINUE
      WRITE( NOUT, FMT = 9985 )
      GO TO 230

  220 CONTINUE
      WRITE( NOUT, FMT = 9991 )

  230 CONTINUE
      IF( TRACE ) CLOSE ( NTRA )
      CLOSE ( NOUT )
      STOP

10002 FORMAT( ' COLUMN-MAJOR AND ROW-MAJOR DATA LAYOUTS ARE TESTED' )
10001 FORMAT(' ROW-MAJOR DATA LAYOUT IS TESTED' )
10000 FORMAT(' COLUMN-MAJOR DATA LAYOUT IS TESTED' )
 9999 FORMAT(' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     $      'S THAN', F8.2 )
 9998 FORMAT(' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 )
 9997 FORMAT(' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ',
     $      'THAN ', I2 )
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 )
 9995 FORMAT(' TESTS OF THE COMPLEX          LEVEL 3 BLAS', //' THE F',
     $      'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9994 FORMAT( '   FOR N              ', 9I6 )
 9993 FORMAT( '   FOR ALPHA          ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9992 FORMAT( '   FOR BETA           ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9991 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM',
     $      /' ******* TESTS ABANDONED *******' )
 9990 FORMAT(' SUBPROGRAM NAME ', A12,' NOT RECOGNIZED', /' ******* T',
     $      'ESTS ABANDONED *******' )
 9989 FORMAT(' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     $      'ATED WRONGLY.', /' CMMCH WAS CALLED WITH TRANSA = ', A1,
     $      'AND TRANSB = ', A1, /' AND RETURNED SAME = ', L1, ' AND ',
     $    ' ERR = ', F12.3, '.', /' THIS MAY BE DUE TO FAULTS IN THE ',
     $     'ARITHMETIC OR THE COMPILER.', /' ******* TESTS ABANDONED ',
     $      '*******' )
 9988 FORMAT( A12,L2 )
 9987 FORMAT( 1X, A12,' WAS NOT TESTED' )
 9986 FORMAT( /' END OF TESTS' )
 9985 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' )
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )

      // End of CBLAT3.

      }
      SUBROUTINE CCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER )

*  Tests CGEMM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS
      REAL               ERR, ERRMAX
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'NTC'/
      // .. Executable Statements ..

      NARGS = 13
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 110 IM = 1, NIDIM
         M = IDIM( IM )

         DO 100 IN = 1, NIDIM
            N = IDIM( IN )
            // Set LDC to 1 more than minimum value if room.
            LDC = M
            IF( LDC.LT.NMAX ) LDC = LDC + 1
            // Skip tests if not enough room.
            IF( LDC.GT.NMAX ) GO TO 100
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0

            DO 90 IK = 1, NIDIM
               K = IDIM( IK )

               DO 80 ICA = 1, 3
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
                  IF( LDA.LT.NMAX ) LDA = LDA + 1
                  // Skip tests if not enough room.
                  IF( LDA.GT.NMAX ) GO TO 80
                  LAA = LDA*NA

                  // Generate the matrix A.

                  CALL CMAKE( 'ge', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO )

                  DO 70 ICB = 1, 3
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
                     IF( LDB.LT.NMAX ) LDB = LDB + 1
                     // Skip tests if not enough room.
                     IF( LDB.GT.NMAX ) GO TO 70
                     LBB = LDB*NB

                     // Generate the matrix B.

                     CALL CMAKE( 'ge', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO )

                     DO 60 IA = 1, NALF
                        ALPHA = ALF( IA )

                        DO 50 IB = 1, NBET
                           BETA = BET( IB )

                           // Generate the matrix C.

                           CALL CMAKE( 'ge', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO )

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           TRANAS = TRANSA
                           TRANBS = TRANSB
                           MS = M
                           NS = N
                           KS = K
                           ALS = ALPHA
                           DO 10 I = 1, LAA
                              AS( I ) = AA( I )
   10                      CONTINUE
                           LDAS = LDA
                           DO 20 I = 1, LBB
                              BS( I ) = BB( I )
   20                      CONTINUE
                           LDBS = LDB
                           BLS = BETA
                           DO 30 I = 1, LCC
                              CS( I ) = CC( I )
   30                      CONTINUE
                           LDCS = LDC

                           // Call the subroutine.

                           IF( TRACE ) CALL CPRCN1(NTRA, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC)
                           IF( REWI ) REWIND NTRA                            CALL CCGEMM( IORDER, TRANSA, TRANSB, M, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC )

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
                           ISAME( 7 ) = LCE( AS, AA, LAA )
                           ISAME( 8 ) = LDAS.EQ.LDA
                           ISAME( 9 ) = LCE( BS, BB, LBB )
                           ISAME( 10 ) = LDBS.EQ.LDB
                           ISAME( 11 ) = BLS.EQ.BETA
                           if ( NULL ) {
                              ISAME( 12 ) = LCE( CS, CC, LCC )
                           } else {
                             ISAME( 12 ) = LCERES( 'ge', ' ', M, N, CS, CC, LDC )
                           }
                           ISAME( 13 ) = LDCS.EQ.LDC

                           // If data was incorrectly changed, report
                           // and return.

                           SAME = .TRUE.
                           DO 40 I = 1, NARGS
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                      CONTINUE
                           if ( .NOT.SAME ) {
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           if ( .NOT.NULL ) {

                              // Check the result.

                             CALL CMMCH( TRANSA, TRANSB, M, N, K, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 120
                           }

   50                   CONTINUE

   60                CONTINUE

   70             CONTINUE

   80          CONTINUE

   90       CONTINUE

  100    CONTINUE

  110 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10000 )SNAME, NC
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10001 )SNAME, NC
      } else {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX
      }
      GO TO 130

  120 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      CALL CPRCN1(NOUT, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC)

  130 CONTINUE
      RETURN

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A12,'(''', A1, ''',''', A1, ''',',
     $     3( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3,
     $     ',(', F4.1, ',', F4.1, '), C,', I3, ').' )
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )

      // End of CCHK1.

      }

      SUBROUTINE CPRCN1(NOUT, NC, SNAME, IORDER, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC)
      int              NOUT, NC, IORDER, M, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA
      String           TRANSA, TRANSB;
      String           SNAME;
      String           CRC, CTA,CTB;

      if (TRANSA.EQ.'N') {
         CTA = '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CTA = '    CblasTrans'
      } else {
         CTA = 'CblasConjTrans'
      }
      if (TRANSB.EQ.'N') {
         CTB = '  CblasNoTrans'
      } else if (TRANSB.EQ.'T') {
         CTB = '    CblasTrans'
      } else {
         CTB = 'CblasConjTrans'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CTA,CTB
      WRITE(NOUT, FMT = 9994)M, N, K, ALPHA, LDA, LDB, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',')
 9994 FORMAT( 10X, 3( I3, ',' ) ,' (', F4.1,',',F4.1,') , A,',
     $ I3, ', B,', I3, ', (', F4.1,',',F4.1,') , C,', I3, ').' )
      }

      SUBROUTINE CCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER )

*  Tests CHEMM and CSYMM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS
      REAL               ERR, ERRMAX
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHS/'LR'/, ICHU/'UL'/
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ).EQ.'he'

      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 100 IM = 1, NIDIM
         M = IDIM( IM )

         DO 90 IN = 1, NIDIM
            N = IDIM( IN )
            // Set LDC to 1 more than minimum value if room.
            LDC = M
            IF( LDC.LT.NMAX ) LDC = LDC + 1
            // Skip tests if not enough room.
            IF( LDC.GT.NMAX ) GO TO 90
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0
            // Set LDB to 1 more than minimum value if room.
            LDB = M
            IF( LDB.LT.NMAX ) LDB = LDB + 1
            // Skip tests if not enough room.
            IF( LDB.GT.NMAX ) GO TO 90
            LBB = LDB*N

            // Generate the matrix B.

            CALL CMAKE( 'ge', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO )

            DO 80 ICS = 1, 2
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'

               if ( LEFT ) {
                  NA = M
               } else {
                  NA = N
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 80
               LAA = LDA*NA

               DO 70 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )

                  // Generate the hermitian or symmetric matrix A.

                  CALL CMAKE(SNAME( 8: 9 ), UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO )

                  DO 60 IA = 1, NALF
                     ALPHA = ALF( IA )

                     DO 50 IB = 1, NBET
                        BETA = BET( IB )

                        // Generate the matrix C.

                        CALL CMAKE( 'ge', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO )

                        NC = NC + 1

                        // Save every datum before calling the
                        // subroutine.

                        SIDES = SIDE
                        UPLOS = UPLO
                        MS = M
                        NS = N
                        ALS = ALPHA
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        DO 20 I = 1, LBB
                           BS( I ) = BB( I )
   20                   CONTINUE
                        LDBS = LDB
                        BLS = BETA
                        DO 30 I = 1, LCC
                           CS( I ) = CC( I )
   30                   CONTINUE
                        LDCS = LDC

                        // Call the subroutine.

                        IF( TRACE ) CALL CPRCN2(NTRA, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC)
                        IF( REWI ) REWIND NTRA
                        if ( CONJ ) {
                           CALL CCHEMM( IORDER, SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC )
                        } else {
                           CALL CCSYMM( IORDER, SIDE, UPLO, M, N, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC )
                        }

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
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LCE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        ISAME( 10 ) = BLS.EQ.BETA
                        if ( NULL ) {
                           ISAME( 11 ) = LCE( CS, CC, LCC )
                        } else {
                           ISAME( 11 ) = LCERES( 'ge', ' ', M, N, CS, CC, LDC )
                        }
                        ISAME( 12 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        DO 40 I = 1, NARGS
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                   CONTINUE
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 110
                        }

                        if ( .NOT.NULL ) {

                           // Check the result.

                           if ( LEFT ) {
                              CALL CMMCH( 'N', 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                           } else {
                              CALL CMMCH( 'N', 'N', M, N, N, ALPHA, B, NMAX, A, NMAX, BETA, C, NMAX, CT, G, CC, LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and
                           // return.
                           IF( FATAL ) GO TO 110
                        }

   50                CONTINUE

   60             CONTINUE

   70          CONTINUE

   80       CONTINUE

   90    CONTINUE

  100 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10000 )SNAME, NC
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10001 )SNAME, NC
      } else {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX
      }
      GO TO 120

  110 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      CALL CPRCN2(NOUT, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC)

  120 CONTINUE
      RETURN

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' )
 9995 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1,
     $      ',', F4.1, '), C,', I3, ')    .' )
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )

      // End of CCHK2.

      }

      SUBROUTINE CPRCN2(NOUT, NC, SNAME, IORDER, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC)
      int              NOUT, NC, IORDER, M, N, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA
      String           SIDE, UPLO;
      String           SNAME;
      String           CRC, CS,CU;

      if (SIDE.EQ.'L') {
         CS =  '     CblasLeft'
      } else {
         CS =  '    CblasRight'
      }
      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CS,CU
      WRITE(NOUT, FMT = 9994)M, N, ALPHA, LDA, LDB, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',')
 9994 FORMAT( 10X, 2( I3, ',' ),' (',F4.1,',',F4.1, '), A,', I3,
     $ ', B,', I3, ', (',F4.1,',',F4.1, '), ', 'C,', I3, ').' )
      }

      SUBROUTINE CCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NMAX, A, AA, AS, B, BB, BS, CT, G, C, IORDER )

*  Tests CTRMM and CTRSM.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RZERO
      const              RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                NALF, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CT( NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS
      REAL               ERR, ERRMAX
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA              ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/, ICHS/'LR'/
      // .. Executable Statements ..

      NARGS = 11
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
      // Set up zero matrix for CMMCH.
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            C( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE

      DO 140 IM = 1, NIDIM
         M = IDIM( IM )

         DO 130 IN = 1, NIDIM
            N = IDIM( IN )
            // Set LDB to 1 more than minimum value if room.
            LDB = M
            IF( LDB.LT.NMAX ) LDB = LDB + 1
            // Skip tests if not enough room.
            IF( LDB.GT.NMAX ) GO TO 130
            LBB = LDB*N
            NULL = M.LE.0.OR.N.LE.0

            DO 120 ICS = 1, 2
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'
               if ( LEFT ) {
                  NA = M
               } else {
                  NA = N
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = NA
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 130
               LAA = LDA*NA

               DO 110 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )

                  DO 100 ICT = 1, 3
                     TRANSA = ICHT( ICT: ICT )

                     DO 90 ICD = 1, 2
                        DIAG = ICHD( ICD: ICD )

                        DO 80 IA = 1, NALF
                           ALPHA = ALF( IA )

                           // Generate the matrix A.

                           CALL CMAKE( 'tr', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO )

                           // Generate the matrix B.

                           CALL CMAKE( 'ge', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO )

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
                           DO 30 I = 1, LAA
                              AS( I ) = AA( I )
   30                      CONTINUE
                           LDAS = LDA
                           DO 40 I = 1, LBB
                              BS( I ) = BB( I )
   40                      CONTINUE
                           LDBS = LDB

                           // Call the subroutine.

                           if ( SNAME( 10: 11 ).EQ.'mm' ) {
                              IF( TRACE ) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB)
                              IF( REWI ) REWIND NTRA                               CALL CCTRMM(IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB )
                           } else if ( SNAME( 10: 11 ).EQ.'sm' ) {
                              IF( TRACE ) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB)
                              IF( REWI ) REWIND NTRA                               CALL CCTRSM(IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA, LDA, BB, LDB )
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
                           ISAME( 8 ) = LCE( AS, AA, LAA )
                           ISAME( 9 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 10 ) = LCE( BS, BB, LBB )
                           } else {
                             ISAME( 10 ) = LCERES( 'ge', ' ', M, N, BS, BB, LDB )
                           }
                           ISAME( 11 ) = LDBS.EQ.LDB

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = .TRUE.
                           DO 50 I = 1, NARGS
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   50                      CONTINUE
                           if ( .NOT.SAME ) {
                              FATAL = .TRUE.
                              GO TO 150
                           }

                           if ( .NOT.NULL ) {
                              if ( SNAME( 10: 11 ).EQ.'mm' ) {

                                 // Check the result.

                                 if ( LEFT ) {
                                   CALL CMMCH( TRANSA, 'N', M, N, M, ALPHA, A, NMAX, B, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .TRUE. )
                                 } else {
                                    CALL CMMCH( 'N', TRANSA, M, N, N, ALPHA, B, NMAX, A, NMAX, ZERO, C, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .TRUE. )
                                 }
                              } else if ( SNAME( 10: 11 ).EQ.'sm' ) {

                                 // Compute approximation to original
                                 // matrix.

                                 DO 70 J = 1, N
                                    DO 60 I = 1, M
                                       C( I, J ) = BB( I + ( J - 1 )* LDB )                                        BB( I + ( J - 1 )*LDB ) = ALPHA* B( I, J )
   60                               CONTINUE
   70                            CONTINUE

                                 if ( LEFT ) {
                                    CALL CMMCH( TRANSA, 'N', M, N, M, ONE, A, NMAX, C, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .FALSE. )
                                 } else {
                                    CALL CMMCH( 'N', TRANSA, M, N, N, ONE, C, NMAX, A, NMAX, ZERO, B, NMAX, CT, G, BB, LDB, EPS, ERR, FATAL, NOUT, .FALSE. )
                                 }
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 150
                           }

   80                   CONTINUE

   90                CONTINUE

  100             CONTINUE

  110          CONTINUE

  120       CONTINUE

  130    CONTINUE

  140 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10000 )SNAME, NC
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10001 )SNAME, NC
      } else {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX
      }
      GO TO 160

  150 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      IF( TRACE ) CALL CPRCN3( NTRA, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB)

  160 CONTINUE
      RETURN

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9996 FORMAT(' ******* ', A12,' FAILED ON CALL NUMBER:' )
 9995 FORMAT(1X, I6, ': ', A12,'(', 4( '''', A1, ''',' ), 2( I3, ',' ),
     $     '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ')         ',
     $      '      .' )
 9994 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )

      // End of CCHK3.

      }

      SUBROUTINE CPRCN3(NOUT, NC, SNAME, IORDER, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB)
      int              NOUT, NC, IORDER, M, N, LDA, LDB;
      COMPLEX          ALPHA
      String           SIDE, UPLO, TRANSA, DIAG;
      String           SNAME;
      String           CRC, CS, CU, CA, CD;

      if (SIDE.EQ.'L') {
         CS =  '     CblasLeft'
      } else {
         CS =  '    CblasRight'
      }
      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (TRANSA.EQ.'N') {
         CA =  '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CA =  '    CblasTrans'
      } else {
         CA =  'CblasConjTrans'
      }
      if (DIAG.EQ.'N') {
         CD =  '  CblasNonUnit'
      } else {
         CD =  '     CblasUnit'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC,SNAME,CRC, CS,CU
      WRITE(NOUT, FMT = 9994)CA, CD, M, N, ALPHA, LDA, LDB

 9995 FORMAT( 1X, I6, ': ', A12,'(', A14, ',', A14, ',', A14, ',')
 9994 FORMAT( 10X, 2( A14, ',') , 2( I3, ',' ), ' (', F4.1, ',',
     $    F4.1, '), A,', I3, ', B,', I3, ').' )
      }

      SUBROUTINE CCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC, CS, CT, G, IORDER )

*  Tests CHERK and CSYRK.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RONE, RZERO
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), B( NMAX, NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS
      REAL               ERR, ERRMAX, RALPHA, RALS, RBETA, RBETS
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NC'/, ICHU/'UL'/
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ).EQ.'he'

      NARGS = 10
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 100 IN = 1, NIDIM
         N = IDIM( IN )
         // Set LDC to 1 more than minimum value if room.
         LDC = N
         IF( LDC.LT.NMAX ) LDC = LDC + 1
         // Skip tests if not enough room.
         IF( LDC.GT.NMAX ) GO TO 100
         LCC = LDC*N

         DO 90 IK = 1, NIDIM
            K = IDIM( IK )

            DO 80 ICT = 1, 2
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'C'
               IF( TRAN.AND..NOT.CONJ ) TRANS = 'T'
               if ( TRAN ) {
                  MA = K
                  NA = N
               } else {
                  MA = N
                  NA = K
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 80
               LAA = LDA*NA

               // Generate the matrix A.

               CALL CMAKE( 'ge', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO )

               DO 70 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'

                  DO 60 IA = 1, NALF
                     ALPHA = ALF( IA )
                     if ( CONJ ) {
                        RALPHA = REAL( ALPHA )
                        ALPHA = CMPLX( RALPHA, RZERO )
                     }

                     DO 50 IB = 1, NBET
                        BETA = BET( IB )
                        if ( CONJ ) {
                           RBETA = REAL( BETA )
                           BETA = CMPLX( RBETA, RZERO )
                        }
                        NULL = N.LE.0
                        IF( CONJ ) NULL = NULL.OR.( ( K.LE.0.OR.RALPHA.EQ. RZERO ).AND.RBETA.EQ.RONE )

                        // Generate the matrix C.

                        CALL CMAKE( SNAME( 8: 9 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO )

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        if ( CONJ ) {
                           RALS = RALPHA
                        } else {
                           ALS = ALPHA
                        }
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        if ( CONJ ) {
                           RBETS = RBETA
                        } else {
                           BETS = BETA
                        }
                        DO 20 I = 1, LCC
                           CS( I ) = CC( I )
   20                   CONTINUE
                        LDCS = LDC

                        // Call the subroutine.

                        if ( CONJ ) {
                           IF( TRACE ) CALL CPRCN6( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC)
                           IF( REWI ) REWIND NTRA                            CALL CCHERK( IORDER, UPLO, TRANS, N, K, RALPHA, AA, LDA, RBETA, CC, LDC )
                        } else {
                           IF( TRACE ) CALL CPRCN4( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC)
                           IF( REWI ) REWIND NTRA                            CALL CCSYRK( IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BETA, CC, LDC )
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9992 )
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS.EQ.UPLO
                        ISAME( 2 ) = TRANSS.EQ.TRANS
                        ISAME( 3 ) = NS.EQ.N
                        ISAME( 4 ) = KS.EQ.K
                        if ( CONJ ) {
                           ISAME( 5 ) = RALS.EQ.RALPHA
                        } else {
                           ISAME( 5 ) = ALS.EQ.ALPHA
                        }
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        if ( CONJ ) {
                           ISAME( 8 ) = RBETS.EQ.RBETA
                        } else {
                           ISAME( 8 ) = BETS.EQ.BETA
                        }
                        if ( NULL ) {
                           ISAME( 9 ) = LCE( CS, CC, LCC )
                        } else {
                           ISAME( 9 ) = LCERES( SNAME( 8: 9 ), UPLO, N, N, CS, CC, LDC )
                        }
                        ISAME( 10 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        DO 30 I = 1, NARGS
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   30                   CONTINUE
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        if ( .NOT.NULL ) {

                           // Check the result column by column.

                           if ( CONJ ) {
                              TRANST = 'C'
                           } else {
                              TRANST = 'T'
                           }
                           JC = 1
                           DO 40 J = 1, N
                              if ( UPPER ) {
                                 JJ = 1
                                 LJ = J
                              } else {
                                 JJ = J
                                 LJ = N - J + 1
                              }
                              if ( TRAN ) {
                                 CALL CMMCH( TRANST, 'N', LJ, 1, K, ALPHA, A( 1, JJ ), NMAX, A( 1, J ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                              } else {
                                 CALL CMMCH( 'N', TRANST, LJ, 1, K, ALPHA, A( JJ, 1 ), NMAX, A( J, 1 ), NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC
                              } else {
                                 JC = JC + LDC + 1
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 110
   40                      CONTINUE
                        }

   50                CONTINUE

   60             CONTINUE

   70          CONTINUE

   80       CONTINUE

   90    CONTINUE

  100 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10000 )SNAME, NC
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10001 )SNAME, NC
      } else {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX
      }
      GO TO 130

  110 CONTINUE
      IF( N.GT.1 ) WRITE( NOUT, FMT = 9995 )J

  120 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( CONJ ) {
      CALL CPRCN6( NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, RALPHA, LDA, rBETA, LDC)
      } else {
      CALL CPRCN4( NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC)
      }

  130 CONTINUE
      RETURN

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $     F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ')               ',
     $      '          .' )
 9993 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, ') , A,', I3, ',(', F4.1, ',', F4.1,
     $      '), C,', I3, ')          .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )

      // End of CCHK4.

      }

      SUBROUTINE CPRCN4(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, BETA, LDC)
      int              NOUT, NC, IORDER, N, K, LDA, LDC;
      COMPLEX          ALPHA, BETA
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (TRANSA.EQ.'N') {
         CA =  '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CA =  '    CblasTrans'
      } else {
         CA =  'CblasConjTrans'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') )
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1 ,'), A,',
     $        I3, ', (', F4.1,',', F4.1, '), C,', I3, ').' )
      }


      SUBROUTINE CPRCN6(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, BETA, LDC)
      int              NOUT, NC, IORDER, N, K, LDA, LDC;
      REAL             ALPHA, BETA
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (TRANSA.EQ.'N') {
         CA =  '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CA =  '    CblasTrans'
      } else {
         CA =  'CblasConjTrans'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') )
 9994 FORMAT( 10X, 2( I3, ',' ),
     $      F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ').' )
      }

      SUBROUTINE CCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W, IORDER )

*  Tests CHER2K and CSYR2K.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      REAL               RONE, RZERO
      const              RONE = 1.0, RZERO = 0.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                NALF, NBET, NIDIM, NMAX, NOUT, NTRA, IORDER;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX            AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ), CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ), W( 2*NMAX )
      REAL               G( NMAX )
      int                IDIM( NIDIM );
      // .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS
      REAL               ERR, ERRMAX, RBETA, RBETS
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
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHT/'NC'/, ICHU/'UL'/
      // .. Executable Statements ..
      CONJ = SNAME( 8: 9 ).EQ.'he'

      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 130 IN = 1, NIDIM
         N = IDIM( IN )
         // Set LDC to 1 more than minimum value if room.
         LDC = N
         IF( LDC.LT.NMAX ) LDC = LDC + 1
         // Skip tests if not enough room.
         IF( LDC.GT.NMAX ) GO TO 130
         LCC = LDC*N

         DO 120 IK = 1, NIDIM
            K = IDIM( IK )

            DO 110 ICT = 1, 2
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'C'
               IF( TRAN.AND..NOT.CONJ ) TRANS = 'T'
               if ( TRAN ) {
                  MA = K
                  NA = N
               } else {
                  MA = N
                  NA = K
               }
               // Set LDA to 1 more than minimum value if room.
               LDA = MA
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 110
               LAA = LDA*NA

               // Generate the matrix A.

               if ( TRAN ) {
                  CALL CMAKE( 'ge', ' ', ' ', MA, NA, AB, 2*NMAX, AA, LDA, RESET, ZERO )
               } else {
                 CALL CMAKE( 'ge', ' ', ' ', MA, NA, AB, NMAX, AA, LDA, RESET, ZERO )
               }

               // Generate the matrix B.

               LDB = LDA
               LBB = LAA
               if ( TRAN ) {
                  CALL CMAKE( 'ge', ' ', ' ', MA, NA, AB( K + 1 ), 2*NMAX, BB, LDB, RESET, ZERO )
               } else {
                  CALL CMAKE( 'ge', ' ', ' ', MA, NA, AB( K*NMAX + 1 ), NMAX, BB, LDB, RESET, ZERO )
               }

               DO 100 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'

                  DO 90 IA = 1, NALF
                     ALPHA = ALF( IA )

                     DO 80 IB = 1, NBET
                        BETA = BET( IB )
                        if ( CONJ ) {
                           RBETA = REAL( BETA )
                           BETA = CMPLX( RBETA, RZERO )
                        }
                        NULL = N.LE.0
                        IF( CONJ ) NULL = NULL.OR.( ( K.LE.0.OR.ALPHA.EQ. ZERO ).AND.RBETA.EQ.RONE )

                        // Generate the matrix C.

                        CALL CMAKE( SNAME( 8: 9 ), UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO )

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        ALS = ALPHA
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        DO 20 I = 1, LBB
                           BS( I ) = BB( I )
   20                   CONTINUE
                        LDBS = LDB
                        if ( CONJ ) {
                           RBETS = RBETA
                        } else {
                           BETS = BETA
                        }
                        DO 30 I = 1, LCC
                           CS( I ) = CC( I )
   30                   CONTINUE
                        LDCS = LDC

                        // Call the subroutine.

                        if ( CONJ ) {
                           IF( TRACE ) CALL CPRCN7( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC)
                           IF( REWI ) REWIND NTRA                            CALL CCHER2K( IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, RBETA, CC, LDC )
                        } else {
                           IF( TRACE ) CALL CPRCN5( NTRA, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC)
                           IF( REWI ) REWIND NTRA                            CALL CCSYR2K( IORDER, UPLO, TRANS, N, K, ALPHA, AA, LDA, BB, LDB, BETA, CC, LDC )
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9992 )
                           FATAL = .TRUE.
                           GO TO 150
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLOS.EQ.UPLO
                        ISAME( 2 ) = TRANSS.EQ.TRANS
                        ISAME( 3 ) = NS.EQ.N
                        ISAME( 4 ) = KS.EQ.K
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LCE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        if ( CONJ ) {
                           ISAME( 10 ) = RBETS.EQ.RBETA
                        } else {
                           ISAME( 10 ) = BETS.EQ.BETA
                        }
                        if ( NULL ) {
                           ISAME( 11 ) = LCE( CS, CC, LCC )
                        } else {
                           ISAME( 11 ) = LCERES( 'he', UPLO, N, N, CS, CC, LDC )
                        }
                        ISAME( 12 ) = LDCS.EQ.LDC

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        DO 40 I = 1, NARGS
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                   CONTINUE
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 150
                        }

                        if ( .NOT.NULL ) {

                           // Check the result column by column.

                           if ( CONJ ) {
                              TRANST = 'C'
                           } else {
                              TRANST = 'T'
                           }
                           JJAB = 1
                           JC = 1
                           DO 70 J = 1, N
                              if ( UPPER ) {
                                 JJ = 1
                                 LJ = J
                              } else {
                                 JJ = J
                                 LJ = N - J + 1
                              }
                              if ( TRAN ) {
                                 DO 50 I = 1, K
                                    W( I ) = ALPHA*AB( ( J - 1 )*2* NMAX + K + I )
                                    if ( CONJ ) {
                                       W( K + I ) = CONJG( ALPHA )* AB( ( J - 1 )*2* NMAX + I )
                                    } else {
                                       W( K + I ) = ALPHA* AB( ( J - 1 )*2* NMAX + I )
                                    }
   50                            CONTINUE
                                 CALL CMMCH( TRANST, 'N', LJ, 1, 2*K, ONE, AB( JJAB ), 2*NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                              } else {
                                 DO 60 I = 1, K
                                    if ( CONJ ) {
                                       W( I ) = ALPHA*CONJG( AB( ( K + I - 1 )*NMAX + J ) )                                        W( K + I ) = CONJG( ALPHA* AB( ( I - 1 )*NMAX + J ) )
                                    } else {
                                       W( I ) = ALPHA*AB( ( K + I - 1 )* NMAX + J )                                        W( K + I ) = ALPHA* AB( ( I - 1 )*NMAX + J )
                                    }
   60                            CONTINUE
                                 CALL CMMCH( 'N', 'N', LJ, 1, 2*K, ONE, AB( JJ ), NMAX, W, 2*NMAX, BETA, C( JJ, J ), NMAX, CT, G, CC( JC ), LDC, EPS, ERR, FATAL, NOUT, .TRUE. )
                              }
                              if ( UPPER ) {
                                 JC = JC + LDC
                              } else {
                                 JC = JC + LDC + 1
                                 IF( TRAN ) JJAB = JJAB + 2*NMAX
                              }
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 140
   70                      CONTINUE
                        }

   80                CONTINUE

   90             CONTINUE

  100          CONTINUE

  110       CONTINUE

  120    CONTINUE

  130 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10000 )SNAME, NC
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10001 )SNAME, NC
      } else {
         IF ( IORDER.EQ.0) WRITE( NOUT, FMT = 10002 )SNAME, NC, ERRMAX
         IF ( IORDER.EQ.1) WRITE( NOUT, FMT = 10003 )SNAME, NC, ERRMAX
      }
      GO TO 160

  140 CONTINUE
      IF( N.GT.1 ) WRITE( NOUT, FMT = 9995 )J

  150 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( CONJ ) {
         CALL CPRCN7( NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC)
      } else {
         CALL CPRCN5( NOUT, NC, SNAME, IORDER, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC)
      }

  160 CONTINUE
      RETURN

10003 FORMAT( ' ', A12,' COMPLETED THE ROW-MAJOR    COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10002 FORMAT( ' ', A12,' COMPLETED THE COLUMN-MAJOR COMPUTATIONAL ',
     $ 'TESTS (', I6, ' CALLS)', /' ******* BUT WITH MAXIMUM TEST ',
     $ 'RATIO ', F8.2, ' - SUSPECT *******' )
10001 FORMAT( ' ', A12,' PASSED THE ROW-MAJOR    COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
10000 FORMAT( ' ', A12,' PASSED THE COLUMN-MAJOR COMPUTATIONAL TESTS',
     $ ' (', I6, ' CALL', 'S)' )
 9998 FORMAT(' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9996 FORMAT( ' ******* ', A12,' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',', F4.1,
     $      ', C,', I3, ')           .' )
 9993 FORMAT(1X, I6, ': ', A12,'(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1,
     $      ',', F4.1, '), C,', I3, ')    .' )
 9992 FORMAT(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )

      // End of CCHK5.

      }

      SUBROUTINE CPRCN5(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, LDB, BETA, LDC)
      int              NOUT, NC, IORDER, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA, BETA
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (TRANSA.EQ.'N') {
         CA =  '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CA =  '    CblasTrans'
      } else {
         CA =  'CblasConjTrans'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, LDB, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') )
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1, '), A,',
     $  I3, ', B', I3, ', (', F4.1, ',', F4.1, '), C,', I3, ').' )
      }


      SUBROUTINE CPRCN7(NOUT, NC, SNAME, IORDER, UPLO, TRANSA, N, K, ALPHA, LDA, LDB, BETA, LDC)
      int              NOUT, NC, IORDER, N, K, LDA, LDB, LDC;
      COMPLEX          ALPHA
      REAL             BETA
      String           UPLO, TRANSA;
      String           SNAME;
      String           CRC, CU, CA;

      if (UPLO.EQ.'U') {
         CU =  '    CblasUpper'
      } else {
         CU =  '    CblasLower'
      }
      if (TRANSA.EQ.'N') {
         CA =  '  CblasNoTrans'
      } else if (TRANSA.EQ.'T') {
         CA =  '    CblasTrans'
      } else {
         CA =  'CblasConjTrans'
      }
      if (IORDER.EQ.1) {
         CRC = ' CblasRowMajor'
      } else {
         CRC = ' CblasColMajor'
      }
      WRITE(NOUT, FMT = 9995)NC, SNAME, CRC, CU, CA
      WRITE(NOUT, FMT = 9994)N, K, ALPHA, LDA, LDB, BETA, LDC

 9995 FORMAT( 1X, I6, ': ', A12,'(', 3( A14, ',') )
 9994 FORMAT( 10X, 2( I3, ',' ), ' (', F4.1, ',', F4.1, '), A,',
     $      I3, ', B', I3, ',', F4.1, ', C,', I3, ').' )
      }

      SUBROUTINE CMAKE(TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, RESET, TRANSL )

*  Generates values for an M by N matrix A.
*  Stores the values in the array AA in the data structure required
*  by the routine, with unwanted elements set to rogue value.

*  TYPE is 'ge', 'he', 'sy' or 'tr'.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      COMPLEX            ROGUE
      const              ROGUE = ( -1.0E10, 1.0E10 ) ;
      REAL               RZERO
      const              RZERO = 0.0 ;
      REAL               RROGUE
      const              RROGUE = -1.0E10 ;
      // .. Scalar Arguments ..
      COMPLEX            TRANSL
      int                LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX            A( NMAX, * ), AA( * )
      // .. Local Scalars ..
      int                I, IBEG, IEND, J, JJ;
      bool               GEN, HER, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      COMPLEX            CBEG
      // EXTERNAL CBEG
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, REAL
      // .. Executable Statements ..
      GEN = TYPE.EQ.'ge'
      HER = TYPE.EQ.'he'
      SYM = TYPE.EQ.'sy'
      TRI = TYPE.EQ.'tr'
      UPPER = ( HER.OR.SYM.OR.TRI ).AND.UPLO.EQ.'U'
      LOWER = ( HER.OR.SYM.OR.TRI ).AND.UPLO.EQ.'L'
      UNIT = TRI.AND.DIAG.EQ.'U'

      // Generate data in array A.

      DO 20 J = 1, N
         DO 10 I = 1, M
            if ( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) ) {
               A( I, J ) = CBEG( RESET ) + TRANSL
               if ( I.NE.J ) {
                  // Set some elements to zero
                  IF( N.GT.3.AND.J.EQ.N/2 ) A( I, J ) = ZERO
                  if ( HER ) {
                     A( J, I ) = CONJG( A( I, J ) )
                  } else if ( SYM ) {
                     A( J, I ) = A( I, J )
                  } else if ( TRI ) {
                     A( J, I ) = ZERO
                  }
               }
            }
   10    CONTINUE
         IF( HER ) A( J, J ) = CMPLX( REAL( A( J, J ) ), RZERO )          IF( TRI ) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE
   20 CONTINUE

      // Store elements in array AS in data structure required by routine.

      if ( TYPE.EQ.'ge' ) {
         DO 50 J = 1, N
            DO 30 I = 1, M
               AA( I + ( J - 1 )*LDA ) = A( I, J )
   30       CONTINUE
            DO 40 I = M + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
   40       CONTINUE
   50    CONTINUE
      } else if ( TYPE.EQ.'he'.OR.TYPE.EQ.'sy'.OR.TYPE.EQ.'tr' ) {
         DO 90 J = 1, N
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
            DO 60 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
   60       CONTINUE
            DO 70 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I, J )
   70       CONTINUE
            DO 80 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
   80       CONTINUE
            if ( HER ) {
               JJ = J + ( J - 1 )*LDA
               AA( JJ ) = CMPLX( REAL( AA( JJ ) ), RROGUE )
            }
   90    CONTINUE
      }
      RETURN

      // End of CMAKE.

      }
      SUBROUTINE CMMCH(TRANSA, TRANSB, M, N, KK, ALPHA, A, LDA, B, LDB, BETA, C, LDC, CT, G, CC, LDCC, EPS, ERR, FATAL, NOUT, MV )

*  Checks the results of the computational tests.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO, RONE
      const              RZERO = 0.0, RONE = 1.0 ;
      // .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      REAL               EPS, ERR
      int                KK, LDA, LDB, LDC, LDCC, M, N, NOUT;
      bool               FATAL, MV;
      String             TRANSA, TRANSB;
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), CC( LDCC, * ), CT( * )
      REAL               G( * )
      // .. Local Scalars ..
      COMPLEX            CL
      REAL               ERRI
      int                I, J, K;
      bool               CTRANA, CTRANB, TRANA, TRANB;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
      // .. Statement Functions ..
      REAL               ABS1
      // .. Statement Function definitions ..
      ABS1( CL ) = ABS( REAL( CL ) ) + ABS( AIMAG( CL ) )
      // .. Executable Statements ..
      TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'
      TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'
      CTRANA = TRANSA.EQ.'C'
      CTRANB = TRANSB.EQ.'C'

      // Compute expected result, one column at a time, in CT using data
      // in A, B and C.
      // Compute gauges in G.

      DO 220 J = 1, N

         DO 10 I = 1, M
            CT( I ) = ZERO
            G( I ) = RZERO
   10    CONTINUE
         if ( .NOT.TRANA.AND..NOT.TRANB ) {
            DO 30 K = 1, KK
               DO 20 I = 1, M
                  CT( I ) = CT( I ) + A( I, K )*B( K, J )
                  G( I ) = G( I ) + ABS1( A( I, K ) )*ABS1( B( K, J ) )
   20          CONTINUE
   30       CONTINUE
         } else if ( TRANA.AND..NOT.TRANB ) {
            if ( CTRANA ) {
               DO 50 K = 1, KK
                  DO 40 I = 1, M
                     CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( K, J )
                     G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) )
   40             CONTINUE
   50          CONTINUE
            } else {
               DO 70 K = 1, KK
                  DO 60 I = 1, M
                     CT( I ) = CT( I ) + A( K, I )*B( K, J )
                     G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( K, J ) )
   60             CONTINUE
   70          CONTINUE
            }
         } else if ( .NOT.TRANA.AND.TRANB ) {
            if ( CTRANB ) {
               DO 90 K = 1, KK
                  DO 80 I = 1, M
                     CT( I ) = CT( I ) + A( I, K )*CONJG( B( J, K ) )
                     G( I ) = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) )
   80             CONTINUE
   90          CONTINUE
            } else {
               DO 110 K = 1, KK
                  DO 100 I = 1, M
                     CT( I ) = CT( I ) + A( I, K )*B( J, K )
                     G( I ) = G( I ) + ABS1( A( I, K ) )* ABS1( B( J, K ) )
  100             CONTINUE
  110          CONTINUE
            }
         } else if ( TRANA.AND.TRANB ) {
            if ( CTRANA ) {
               if ( CTRANB ) {
                  DO 130 K = 1, KK
                     DO 120 I = 1, M
                        CT( I ) = CT( I ) + CONJG( A( K, I ) )* CONJG( B( J, K ) )                         G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) )
  120                CONTINUE
  130             CONTINUE
               } else {
                  DO 150 K = 1, KK
                     DO 140 I = 1, M
                       CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( J, K )
                       G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) )
  140                CONTINUE
  150             CONTINUE
               }
            } else {
               if ( CTRANB ) {
                  DO 170 K = 1, KK
                     DO 160 I = 1, M
                       CT( I ) = CT( I ) + A( K, I )*CONJG( B( J, K ) )
                       G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) )
  160                CONTINUE
  170             CONTINUE
               } else {
                  DO 190 K = 1, KK
                     DO 180 I = 1, M
                        CT( I ) = CT( I ) + A( K, I )*B( J, K )
                        G( I ) = G( I ) + ABS1( A( K, I ) )* ABS1( B( J, K ) )
  180                CONTINUE
  190             CONTINUE
               }
            }
         }
         DO 200 I = 1, M
            CT( I ) = ALPHA*CT( I ) + BETA*C( I, J )
            G( I ) = ABS1( ALPHA )*G( I ) + ABS1( BETA )*ABS1( C( I, J ) )
  200    CONTINUE

         // Compute the error ratio for this result.

         ERR = ZERO
         DO 210 I = 1, M
            ERRI = ABS1( CT( I ) - CC( I, J ) )/EPS
            IF( G( I ).NE.RZERO ) ERRI = ERRI/G( I )
            ERR = MAX( ERR, ERRI )
            IF( ERR*SQRT( EPS ).GE.RONE ) GO TO 230
  210    CONTINUE

  220 CONTINUE

      // If the loop completes, all results are at least half accurate.
      GO TO 250

      // Report fatal error.

  230 FATAL = .TRUE.
      WRITE( NOUT, FMT = 9999 )
      DO 240 I = 1, M
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, CT( I ), CC( I, J )
         } else {
            WRITE( NOUT, FMT = 9998 )I, CC( I, J ), CT( I )
         }
  240 CONTINUE
      IF( N.GT.1 ) WRITE( NOUT, FMT = 9997 )J

  250 CONTINUE
      RETURN

 9999 FORMAT(' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',
     $     'F ACCURATE *******', /'                       EXPECTED RE',
     $     'SULT                    COMPUTED RESULT' )
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) )
 9997 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )

      // End of CMMCH.

      }
      bool    FUNCTION LCE( RI, RJ, LR );

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
      COMPLEX            RI( * ), RJ( * )
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      DO 10 I = 1, LR
         IF( RI( I ).NE.RJ( I ) ) GO TO 20
   10 CONTINUE
      LCE = .TRUE.
      GO TO 30
   20 CONTINUE
      LCE = .FALSE.
   30 RETURN

      // End of LCE.

      }
      bool    FUNCTION LCERES( TYPE, UPLO, M, N, AA, AS, LDA );

*  Tests if selected elements in two arrays are equal.

*  TYPE is 'ge' or 'he' or 'sy'.

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
      COMPLEX            AA( LDA, * ), AS( LDA, * )
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO.EQ.'U'
      if ( TYPE.EQ.'ge' ) {
         DO 20 J = 1, N
            DO 10 I = M + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   10       CONTINUE
   20    CONTINUE
      } else if ( TYPE.EQ.'he'.OR.TYPE.EQ.'sy' ) {
         DO 50 J = 1, N
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            DO 30 I = 1, IBEG - 1
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   30       CONTINUE
            DO 40 I = IEND + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   40       CONTINUE
   50    CONTINUE
      }

   60 CONTINUE
      LCERES = .TRUE.
      GO TO 80
   70 CONTINUE
      LCERES = .FALSE.
   80 RETURN

      // End of LCERES.

      }
      COMPLEX FUNCTION CBEG( RESET )

*  Generates complex numbers as pairs of random numbers uniformly
*  distributed between -0.5 and 0.5.

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX
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
      CBEG = CMPLX( ( I - 500 )/1001.0, ( J - 500 )/1001.0 )
      RETURN

      // End of CBEG.

      }
      REAL FUNCTION SDIFF( X, Y )

*  Auxiliary routine for test program for Level 3 Blas.

*  -- Written on 8-February-1989.
      // Jack Dongarra, Argonne National Laboratory.
      // Iain Duff, AERE Harwell.
      // Jeremy Du Croz, Numerical Algorithms Group Ltd.
      // Sven Hammarling, Numerical Algorithms Group Ltd.

      // .. Scalar Arguments ..
      REAL               X, Y
      // .. Executable Statements ..
      SDIFF = X - Y
      RETURN

      // End of SDIFF.

      }
