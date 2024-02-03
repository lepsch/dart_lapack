      void cdrvpp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NVAL( * );
      REAL               RWORK( * ), S( * );
      COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, K, K1, KL, KU, LDA, MODE, N, NERRS, NFACT, NFAIL, NIMAT, NPP, NRUN, NT;
      REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SGET06;
      // EXTERNAL LSAME, CLANHP, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CCOPY, CERRVX, CGET04, CLACPY, CLAIPD, CLAQHP, CLARHS, CLASET, CLATB4, CLATMS, CPPEQU, CPPSV, CPPSVX, CPPT01, CPPT02, CPPT05, CPPTRF, CPPTRI
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
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', FACTS = 'F', 'N', 'E', PACKS = 'C', 'R', EQUEDS = 'N', 'Y';
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Complex precision';
      PATH( 2: 3 ) = 'PP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) CALL CERRVX( PATH, NOUT );
      INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 140
         N = NVAL( IN );
         LDA = max( N, 1 );
         NPP = N*( N+1 ) / 2;
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 130

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 130;

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 5;
            if (ZEROT && N < IMAT-2) GO TO 130;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 120
               UPLO = UPLOS( IUPLO );
               PACKIT = PACKS( IUPLO );

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
               RCONDC = ONE / CNDNUM;

               SRNAMT = 'CLATMS';
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 120;
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1;
                  } else if ( IMAT == 4 ) {
                     IZERO = N;
                  } else {
                     IZERO = N / 2 + 1;
                  }

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO == 1 ) {
                     IOFF = ( IZERO-1 )*IZERO / 2;
                     for (I = 1; I <= IZERO - 1; I++) { // 20
                        A( IOFF+I ) = ZERO;
                     } // 20
                     IOFF = IOFF + IZERO;
                     for (I = IZERO; I <= N; I++) { // 30
                        A( IOFF ) = ZERO;
                        IOFF = IOFF + I;
                     } // 30
                  } else {
                     IOFF = IZERO;
                     for (I = 1; I <= IZERO - 1; I++) { // 40
                        A( IOFF ) = ZERO;
                        IOFF = IOFF + N - I;
                     } // 40
                     IOFF = IOFF - IZERO;
                     for (I = IZERO; I <= N; I++) { // 50
                        A( IOFF+I ) = ZERO;
                     } // 50
                  }
               } else {
                  IZERO = 0;
               }

               // Set the imaginary part of the diagonals.

               if ( IUPLO == 1 ) {
                  claipd(N, A, 2, 1 );
               } else {
                  claipd(N, A, N, -1 );
               }

               // Save a copy of the matrix A in ASAV.

               ccopy(NPP, A, 1, ASAV, 1 );

               for (IEQUED = 1; IEQUED <= 2; IEQUED++) { // 110
                  EQUED = EQUEDS( IEQUED );
                  if ( IEQUED == 1 ) {
                     NFACT = 3;
                  } else {
                     NFACT = 1;
                  }

                  for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 100
                     FACT = FACTS( IFACT );
                     PREFAC = LSAME( FACT, 'F' );
                     NOFACT = LSAME( FACT, 'N' );
                     EQUIL = LSAME( FACT, 'E' );

                     if ( ZEROT ) {
                        if (PREFAC) GO TO 100;
                        RCONDC = ZERO;

                     } else if ( !LSAME( FACT, 'N' ) ) {

                        // Compute the condition number for comparison with
                        // the value returned by CPPSVX (FACT = 'N' reuses
                        // the condition number from the previous iteration
                           // with FACT = 'F').

                        ccopy(NPP, ASAV, 1, AFAC, 1 );
                        if ( EQUIL || IEQUED > 1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           cppequ(UPLO, N, AFAC, S, SCOND, AMAX, INFO );
                           if ( INFO == 0 && N > 0 ) {
                              if (IEQUED > 1) SCOND = ZERO;

                              // Equilibrate the matrix.

                              claqhp(UPLO, N, AFAC, S, SCOND, AMAX, EQUED );
                           }
                        }

                        // Save the condition number of the
                        // non-equilibrated system for use in CGET04.

                        if (EQUIL) ROLDC = RCONDC;

                        // Compute the 1-norm of A.

                        ANORM = CLANHP( '1', UPLO, N, AFAC, RWORK );

                        // Factor the matrix A.

                        cpptrf(UPLO, N, AFAC, INFO );

                        // Form the inverse of A.

                        ccopy(NPP, AFAC, 1, A, 1 );
                        cpptri(UPLO, N, A, INFO );

                        // Compute the 1-norm condition number of A.

                        AINVNM = CLANHP( '1', UPLO, N, A, RWORK );
                        if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                           RCONDC = ONE;
                        } else {
                           RCONDC = ( ONE / ANORM ) / AINVNM;
                        }
                     }

                     // Restore the matrix A.

                     ccopy(NPP, ASAV, 1, A, 1 );

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'CLARHS';
                     clarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     clacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT ) {

                        // --- Test CPPSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        ccopy(NPP, A, 1, AFAC, 1 );
                        clacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'CPPSV ';
                        cppsv(UPLO, N, NRHS, AFAC, X, LDA, INFO );

                        // Check error code from CPPSV .

                        if ( INFO != IZERO ) {
                           alaerh(PATH, 'CPPSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 70;
                        } else if ( INFO != 0 ) {
                           GO TO 70;
                        }

                        // Reconstruct matrix from factors and compute
                        // residual.

                        cppt01(UPLO, N, A, AFAC, RWORK, RESULT( 1 ) );

                        // Compute residual of the computed solution.

                        clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        cppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        NT = 3;

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 60
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                              WRITE( NOUT, FMT = 9999 )'CPPSV ', UPLO, N, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 60
                        NRUN = NRUN + NT;
                        } // 70
                     }

                     // --- Test CPPSVX ---

                     if ( !PREFAC && NPP > 0) CALL CLASET( 'Full', NPP, 1, CMPLX( ZERO ), CMPLX( ZERO ), AFAC, NPP );
                     claset('Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDA );
                     if ( IEQUED > 1 && N > 0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        claqhp(UPLO, N, A, S, SCOND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using CPPSVX.

                     SRNAMT = 'CPPSVX';
                     cppsvx(FACT, UPLO, N, NRHS, A, AFAC, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check the error code from CPPSVX.

                     if ( INFO != IZERO ) {
                        alaerh(PATH, 'CPPSVX', INFO, IZERO, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 90;
                     }

                     if ( INFO == 0 ) {
                        if ( !PREFAC ) {

                           // Reconstruct matrix from factors and compute
                           // residual.

                           cppt01(UPLO, N, A, AFAC, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                           K1 = 1;
                        } else {
                           K1 = 2;
                        }

                        // Compute residual of the computed solution.

                        clacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA );
                        cppt02(UPLO, N, NRHS, ASAV, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        if( NOFACT || ( PREFAC && LSAME( EQUED, 'N' ) ) ) THEN;
                           cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        } else {
                           cget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        cppt05(UPLO, N, NRHS, ASAV, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        K1 = 6;
                     }

                     // Compare RCOND from CPPSVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = SGET06( RCOND, RCONDC );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = K1; K <= 6; K++) { // 80
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'CPPSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT( K );
                           } else {
                              WRITE( NOUT, FMT = 9998 )'CPPSVX', FACT, UPLO, N, IMAT, K, RESULT( K );
                           }
                           NFAIL = NFAIL + 1;
                        }
                     } // 80
                     NRUN = NRUN + 7 - K1;
                     } // 90
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 );
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', type ', I1, ', test(', I1, ')=', G12.5 );
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=', G12.5 );
      return;

      // End of CDRVPP

      }
