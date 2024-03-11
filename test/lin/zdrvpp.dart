      void zdrvpp(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NRHS, final double THRESH, final bool TSTERR, final int NMAX, final Array<Complex> A_, final Array<Complex> AFAC_, final Array<Complex> ASAV_, final Array<Complex> B_, final Array<Complex> BSAV_, final Array<double> X_, final Array<double> XACT_, final int S, final Array<double> WORK_, final Array<double> RWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      bool               DOTYPE( * );
      int                NVAL( * );
      double             RWORK( * ), S( * );
      Complex         A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, K, K1, KL, KU, LDA, MODE, N, NFACT, NIMAT, NPP, NRUN, NT;
      double             AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      String             EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 );
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DGET06, ZLANHP;
      // EXTERNAL lsame, DGET06, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, ZCOPY, ZERRVX, ZGET04, ZLACPY, ZLAIPD, ZLAQHP, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPPEQU, ZPPSV, ZPPSVX, ZPPT01, ZPPT02, ZPPT05, ZPPTRF, ZPPTRI
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', FACTS = 'F', 'N', 'E', PACKS = 'C', 'R', EQUEDS = 'N', 'Y';

      // Initialize constants and the random number seed.

      final PATH = '${'Zomplex precision'[0]}PP';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      // Test the error exits

      if (TSTERR) zerrvx( PATH, NOUT );
      infoc.INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 140
         final N = NVAL[IN];
         final LDA = max( N, 1 );
         NPP = N*( N+1 ) / 2;
         XTYPE = 'N';
            final NIMAT = N <= 0 ? 1 : NTYPES;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 130

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE[IMAT] ) GO TO 130;

            // Skip types 3, 4, or 5 if the matrix size is too small.

            final ZEROT = IMAT >= 3 && IMAT <= 5;
            if (ZEROT && N < IMAT-2) GO TO 130;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 120
               final UPLO = UPLOS[IUPLO - 1];
               PACKIT = PACKS( IUPLO );

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               final (:TYPE,:KL,:KU,:ANORM,:MODE,:CNDNUM,:DIST) = zlatb4(PATH, IMAT, N, N);
               RCONDC = ONE / CNDNUM;

              srnamc.SRNAMT = 'ZLATMS';
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from ZLATMS.

               if ( INFO.value != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 120;
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               final int IZERO;
               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1;
                  } else if ( IMAT == 4 ) {
                     IZERO = N;
                  } else {
                     IZERO = N ~/ 2 + 1;
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
                  final EQUED = EQUEDS[IEQUED - 1];
                  if ( IEQUED == 1 ) {
                     NFACT = 3;
                  } else {
                     NFACT = 1;
                  }

                  for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 100
                     final FACT = FACTS[IFACT - 1];
                     final PREFAC = lsame( FACT, 'F' );
                     final NOFACT = lsame( FACT, 'N' );
                     final EQUIL = lsame( FACT, 'E' );

                     final int IZERO;
                     if ( ZEROT ) {
                        if (PREFAC) GO TO 100;
                        RCONDC = ZERO;

                     } else if ( !lsame( FACT, 'N' ) ) {

                        // Compute the condition number for comparison with
                        // the value returned by ZPPSVX (FACT = 'N' reuses
                        // the condition number from the previous iteration
                        //    with FACT = 'F').

                        zcopy(NPP, ASAV, 1, AFAC, 1 );
                        if ( EQUIL || IEQUED > 1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           zppequ(UPLO, N, AFAC, S, SCOND, AMAX, INFO );
                           if ( INFO.value == 0 && N > 0 ) {
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

                    srnamc.SRNAMT = 'ZLARHS';
                     zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     zlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT ) {

                        // --- Test ZPPSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        zcopy(NPP, A, 1, AFAC, 1 );
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'ZPPSV ';
                        zppsv(UPLO, N, NRHS, AFAC, X, LDA, INFO );

                        // Check error code from ZPPSV .

                        if ( INFO.value != IZERO ) {
                           alaerh(PATH, 'ZPPSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 70;
                        } else if ( INFO.value != 0 ) {
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
                           if ( RESULT[K] >= THRESH ) {
                              if (NFAIL == 0 && NERRS.value == 0) aladhd( NOUT, PATH );
                              NOUT.println( 9999 )'ZPPSV ', UPLO, N, IMAT, K, RESULT[K];
                              NFAIL++;
                           }
                        } // 60
                        NRUN +=  NT;
                        } // 70
                     }

                     // --- Test ZPPSVX ---

                     if ( !PREFAC && NPP > 0) zlaset( 'Full', NPP, 1, Complex.zero, Complex.zero, AFAC, NPP );
                     zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X, LDA );
                     if ( IEQUED > 1 && N > 0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        zlaqhp(UPLO, N, A, S, SCOND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZPPSVX.

                    srnamc.SRNAMT = 'ZPPSVX';
                     zppsvx(FACT, UPLO, N, NRHS, A, AFAC, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check the error code from ZPPSVX.

                     if ( INFO.value != IZERO ) {
                        alaerh(PATH, 'ZPPSVX', INFO, IZERO, FACT + UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 90;
                     }

                     if ( INFO.value == 0 ) {
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

                        if( NOFACT || ( PREFAC && lsame( EQUED, 'N' ) ) ) THEN;
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
                        if ( RESULT[K] >= THRESH ) {
                           if (NFAIL == 0 && NERRS.value == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              NOUT.println( 9997 )'ZPPSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT[K];
                           } else {
                              NOUT.println( 9998 )'ZPPSVX', FACT, UPLO, N, IMAT, K, RESULT[K];
                           }
                           NFAIL++;
                        }
                     } // 80
                     NRUN +=  7 - K1;
                     } // 90
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${}, UPLO=\'${.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${.i1})=${RESULT[].g12_5};
 9998 FORMAT(' ${}, FACT=\'${FACT.a1}\', UPLO=\'${.a1}\', N=${N.i5}, type ${IMAT.i1}, test(${.i1})=${RESULT[].g12_5};
 9997 FORMAT(' ${}, FACT=\'${FACT.a1}\', UPLO=\'${.a1}\', N=${N.i5}, EQUED=\'${EQUED.a1}\', type ${IMAT.i1}, test(${.i1})=${RESULT[].g12_5};
      }
