      void zdrvpp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NVAL( * );
      double             RWORK( * ), S( * );
      Complex         A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
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
      double             AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DGET06, ZLANHP;
      // EXTERNAL LSAME, DGET06, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, ZCOPY, ZERRVX, ZGET04, ZLACPY, ZLAIPD, ZLAQHP, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPPEQU, ZPPSV, ZPPSVX, ZPPT01, ZPPT02, ZPPT05, ZPPTRF, ZPPTRI
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
      // INTRINSIC DCMPLX, MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', FACTS = 'F', 'N', 'E', PACKS = 'C', 'R', EQUEDS = 'N', 'Y';
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'PP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrvx( PATH, NOUT );
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

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
               RCONDC = ONE / CNDNUM;

               SRNAMT = 'ZLATMS';
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from ZLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
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
                        A[IOFF+I] = ZERO;
                     } // 20
                     IOFF = IOFF + IZERO;
                     for (I = IZERO; I <= N; I++) { // 30
                        A[IOFF] = ZERO;
                        IOFF = IOFF + I;
                     } // 30
                  } else {
                     IOFF = IZERO;
                     for (I = 1; I <= IZERO - 1; I++) { // 40
                        A[IOFF] = ZERO;
                        IOFF = IOFF + N - I;
                     } // 40
                     IOFF = IOFF - IZERO;
                     for (I = IZERO; I <= N; I++) { // 50
                        A[IOFF+I] = ZERO;
                     } // 50
                  }
               } else {
                  IZERO = 0;
               }

               // Set the imaginary part of the diagonals.

               if ( IUPLO == 1 ) {
                  zlaipd(N, A, 2, 1 );
               } else {
                  zlaipd(N, A, N, -1 );
               }

               // Save a copy of the matrix A in ASAV.

               zcopy(NPP, A, 1, ASAV, 1 );

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
                        // the value returned by ZPPSVX (FACT = 'N' reuses
                        // the condition number from the previous iteration
                           // with FACT = 'F').

                        zcopy(NPP, ASAV, 1, AFAC, 1 );
                        if ( EQUIL || IEQUED > 1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           zppequ(UPLO, N, AFAC, S, SCOND, AMAX, INFO );
                           if ( INFO == 0 && N > 0 ) {
                              if (IEQUED > 1) SCOND = ZERO;

                              // Equilibrate the matrix.

                              zlaqhp(UPLO, N, AFAC, S, SCOND, AMAX, EQUED );
                           }
                        }

                        // Save the condition number of the
                        // non-equilibrated system for use in ZGET04.

                        if (EQUIL) ROLDC = RCONDC;

                        // Compute the 1-norm of A.

                        ANORM = ZLANHP( '1', UPLO, N, AFAC, RWORK );

                        // Factor the matrix A.

                        zpptrf(UPLO, N, AFAC, INFO );

                        // Form the inverse of A.

                        zcopy(NPP, AFAC, 1, A, 1 );
                        zpptri(UPLO, N, A, INFO );

                        // Compute the 1-norm condition number of A.

                        AINVNM = ZLANHP( '1', UPLO, N, A, RWORK );
                        if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                           RCONDC = ONE;
                        } else {
                           RCONDC = ( ONE / ANORM ) / AINVNM;
                        }
                     }

                     // Restore the matrix A.

                     zcopy(NPP, ASAV, 1, A, 1 );

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'ZLARHS';
                     zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     zlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT ) {

                        // --- Test ZPPSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        zcopy(NPP, A, 1, AFAC, 1 );
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'ZPPSV ';
                        zppsv(UPLO, N, NRHS, AFAC, X, LDA, INFO );

                        // Check error code from ZPPSV .

                        if ( INFO != IZERO ) {
                           alaerh(PATH, 'ZPPSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 70;
                        } else if ( INFO != 0 ) {
                           GO TO 70;
                        }

                        // Reconstruct matrix from factors and compute
                        // residual.

                        zppt01(UPLO, N, A, AFAC, RWORK, RESULT( 1 ) );

                        // Compute residual of the computed solution.

                        zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        zppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        NT = 3;

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 60
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9999 )'ZPPSV ', UPLO, N, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 60
                        NRUN = NRUN + NT;
                        } // 70
                     }

                     // --- Test ZPPSVX ---

                     if ( !PREFAC && NPP > 0) zlaset( 'Full', NPP, 1, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, NPP );
                     zlaset('Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA );
                     if ( IEQUED > 1 && N > 0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        zlaqhp(UPLO, N, A, S, SCOND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZPPSVX.

                     SRNAMT = 'ZPPSVX';
                     zppsvx(FACT, UPLO, N, NRHS, A, AFAC, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check the error code from ZPPSVX.

                     if ( INFO != IZERO ) {
                        alaerh(PATH, 'ZPPSVX', INFO, IZERO, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 90;
                     }

                     if ( INFO == 0 ) {
                        if ( !PREFAC ) {

                           // Reconstruct matrix from factors and compute
                           // residual.

                           zppt01(UPLO, N, A, AFAC, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                           K1 = 1;
                        } else {
                           K1 = 2;
                        }

                        // Compute residual of the computed solution.

                        zlacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA );
                        zppt02(UPLO, N, NRHS, ASAV, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        if( NOFACT || ( PREFAC && LSAME( EQUED, 'N' ) ) ) THEN;
                           zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        } else {
                           zget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        zppt05(UPLO, N, NRHS, ASAV, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        K1 = 6;
                     }

                     // Compare RCOND from ZPPSVX with the computed value
                     // in RCONDC.

                     RESULT[6] = DGET06( RCOND, RCONDC );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = K1; K <= 6; K++) { // 80
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'ZPPSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT( K );
                           } else {
                              WRITE( NOUT, FMT = 9998 )'ZPPSVX', FACT, UPLO, N, IMAT, K, RESULT( K );
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
      }
