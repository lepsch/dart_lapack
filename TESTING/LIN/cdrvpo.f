      SUBROUTINE CDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NVAL( * );
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, K, K1, KL, KU, LDA, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT;
      REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHE, SGET06
      // EXTERNAL LSAME, CLANHE, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CERRVX, CGET04, CLACPY, CLAIPD, CLAQHE, CLARHS, CLASET, CLATB4, CLATMS, CPOEQU, CPOSV, CPOSVX, CPOT01, CPOT02, CPOT05, CPOTRF, CPOTRI, XLAENV
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'Y' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL CERRVX( PATH, NOUT );
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 130
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 120

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            if (ZEROT .AND. N.LT.IMAT-2) GO TO 120;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110
               UPLO = UPLOS( IUPLO )

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'CLATMS'
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 110
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1
                  } else if ( IMAT == 4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO == 1 ) {
                     for (I = 1; I <= IZERO - 1; I++) { // 20
                        A( IOFF+I ) = ZERO
                     } // 20
                     IOFF = IOFF + IZERO
                     for (I = IZERO; I <= N; I++) { // 30
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
                     } // 30
                  } else {
                     IOFF = IZERO
                     for (I = 1; I <= IZERO - 1; I++) { // 40
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
                     } // 40
                     IOFF = IOFF - IZERO
                     for (I = IZERO; I <= N; I++) { // 50
                        A( IOFF+I ) = ZERO
                     } // 50
                  }
               } else {
                  IZERO = 0
               }

               // Set the imaginary part of the diagonals.

               claipd(N, A, LDA+1, 0 );

               // Save a copy of the matrix A in ASAV.

               clacpy(UPLO, N, N, A, LDA, ASAV, LDA );

               for (IEQUED = 1; IEQUED <= 2; IEQUED++) { // 100
                  EQUED = EQUEDS( IEQUED )
                  if ( IEQUED == 1 ) {
                     NFACT = 3
                  } else {
                     NFACT = 1
                  }

                  for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 90
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )

                     if ( ZEROT ) {
                        if (PREFAC) GO TO 90;
                        RCONDC = ZERO

                     } else if ( .NOT.LSAME( FACT, 'N' ) ) {

                        // Compute the condition number for comparison with
                        // the value returned by CPOSVX (FACT = 'N' reuses
                        // the condition number from the previous iteration
                        // with FACT = 'F').

                        clacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );
                        if ( EQUIL .OR. IEQUED.GT.1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           cpoequ(N, AFAC, LDA, S, SCOND, AMAX, INFO );
                           if ( INFO == 0 .AND. N.GT.0 ) {
                              if (IEQUED.GT.1) SCOND = ZERO;

                              // Equilibrate the matrix.

                              claqhe(UPLO, N, AFAC, LDA, S, SCOND, AMAX, EQUED );
                           }
                        }

                        // Save the condition number of the
                        // non-equilibrated system for use in CGET04.

                        if (EQUIL) ROLDC = RCONDC;

                        // Compute the 1-norm of A.

                        ANORM = CLANHE( '1', UPLO, N, AFAC, LDA, RWORK )

                        // Factor the matrix A.

                        cpotrf(UPLO, N, AFAC, LDA, INFO );

                        // Form the inverse of A.

                        clacpy(UPLO, N, N, AFAC, LDA, A, LDA );
                        cpotri(UPLO, N, A, LDA, INFO );

                        // Compute the 1-norm condition number of A.

                        AINVNM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
                        if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                           RCONDC = ONE
                        } else {
                           RCONDC = ( ONE / ANORM ) / AINVNM
                        }
                     }

                     // Restore the matrix A.

                     clacpy(UPLO, N, N, ASAV, LDA, A, LDA );

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'CLARHS'
                     clarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C'
                     clacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT ) {

                        // --- Test CPOSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        clacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                        clacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'CPOSV '
                        cposv(UPLO, N, NRHS, AFAC, LDA, X, LDA, INFO );

                        // Check error code from CPOSV .

                        if ( INFO.NE.IZERO ) {
                           alaerh(PATH, 'CPOSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 70
                        } else if ( INFO.NE.0 ) {
                           GO TO 70
                        }

                        // Reconstruct matrix from factors and compute
                        // residual.

                        cpot01(UPLO, N, A, LDA, AFAC, LDA, RWORK, RESULT( 1 ) );

                        // Compute residual of the computed solution.

                        clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        cpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        NT = 3

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 60
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL == 0 .AND. NERRS == 0) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'CPOSV ', UPLO, N, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1
                           }
                        } // 60
                        NRUN = NRUN + NT
                        } // 70
                     }

                     // --- Test CPOSVX ---

                     if (.NOT.PREFAC) CALL CLASET( UPLO, N, N, CMPLX( ZERO ), CMPLX( ZERO ), AFAC, LDA );
                     claset('Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDA );
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        claqhe(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using CPOSVX.

                     SRNAMT = 'CPOSVX'
                     cposvx(FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check the error code from CPOSVX.

                     if ( INFO.NE.IZERO ) {
                        alaerh(PATH, 'CPOSVX', INFO, IZERO, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 90
                     }

                     if ( INFO == 0 ) {
                        if ( .NOT.PREFAC ) {

                           // Reconstruct matrix from factors and compute
                           // residual.

                           cpot01(UPLO, N, A, LDA, AFAC, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                           K1 = 1
                        } else {
                           K1 = 2
                        }

                        // Compute residual of the computed solution.

                        clacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA );
                        cpot02(UPLO, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN;
                           cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        } else {
                           cget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        cpot05(UPLO, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        K1 = 6
                     }

                     // Compare RCOND from CPOSVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = SGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = K1; K <= 6; K++) { // 80
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL == 0 .AND. NERRS == 0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'CPOSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT( K )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'CPOSVX', FACT, UPLO, N, IMAT, K, RESULT( K )
                           }
                           NFAIL = NFAIL + 1
                        }
                     } // 80
                     NRUN = NRUN + 7 - K1
                  } // 90
               } // 100
            } // 110
         } // 120
      } // 130

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ') =', G12.5 )
      RETURN

      // End of CDRVPO

      }
