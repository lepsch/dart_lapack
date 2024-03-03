      void sdrvge(final int DOTYPE, final int NN, final int NVAL, final int NRHS, final int THRESH, final int TSTERR, final int NMAX, final int A, final int AFAC, final int ASAV, final int B, final int BSAV, final int X, final int XACT, final int S, final Array<double> _WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final int NOUT,) {
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double               THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double               A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), RWORK( * ), S( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      bool               EQUIL, NOFACT, PREFAC, TRFCON, ZEROT;
      String             DIST, EQUED, FACT, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, ITRAN, IZERO, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT;
      double               AINVNM, AMAX, ANORM, ANORMI, ANORMO, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW;
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SGET06, SLAMCH, SLANGE, SLANTR;
      // EXTERNAL lsame, SGET06, SLAMCH, SLANGE, SLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, SERRVX, SGEEQU, SGESV, SGESVX, SGET01, SGET02, SGET04, SGET07, SGETRF, SGETRI, SLACPY, SLAQGE, SLARHS, SLASET, SLATB4, SLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const TRANSS = [ 'N', 'T', 'C' ];
      const FACTS = [ 'F', 'N', 'E' ];
      const EQUEDS = [ 'N', 'R', 'C', 'B' ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'GE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) serrvx( PATH, NOUT );
      INFOT = 0;

      // Set the block size and minimum block size for testing.

      NB = 1;
      NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 90
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 80

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 80;

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT >= 5 && IMAT <= 7;
            if (ZEROT && N < IMAT-4) GO TO 80;

            // Set up parameters with SLATB4 and generate a test matrix
            // with SLATMS.

            slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
            RCONDC = ONE / CNDNUM;

           srnamc.SRNAMT = 'SLATMS';
            slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from SLATMS.

            if ( INFO != 0 ) {
               alaerh(PATH, 'SLATMS', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 80;
            }

            // For types 5-7, zero one or more columns of the matrix to
            // test that INFO is returned correctly.

            if ( ZEROT ) {
               if ( IMAT == 5 ) {
                  IZERO = 1;
               } else if ( IMAT == 6 ) {
                  IZERO = N;
               } else {
                  IZERO = N / 2 + 1;
               }
               IOFF = ( IZERO-1 )*LDA;
               if ( IMAT < 7 ) {
                  for (I = 1; I <= N; I++) { // 20
                     A[IOFF+I] = ZERO;
                  } // 20
               } else {
                  slaset('Full', N, N-IZERO+1, ZERO, ZERO, A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0;
            }

            // Save a copy of the matrix A in ASAV.

            slacpy('Full', N, N, A, LDA, ASAV, LDA );

            for (IEQUED = 1; IEQUED <= 4; IEQUED++) { // 70
               EQUED = EQUEDS( IEQUED );
               if ( IEQUED == 1 ) {
                  NFACT = 3;
               } else {
                  NFACT = 1;
               }

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 60
                  FACT = FACTS( IFACT );
                  PREFAC = lsame( FACT, 'F' );
                  NOFACT = lsame( FACT, 'N' );
                  EQUIL = lsame( FACT, 'E' );

                  if ( ZEROT ) {
                     if (PREFAC) GO TO 60;
                     RCONDO = ZERO;
                     RCONDI = ZERO;

                  } else if ( !NOFACT ) {

                     // Compute the condition number for comparison with
                     // the value returned by SGESVX (FACT = 'N' reuses
                     // the condition number from the previous iteration
                     // with FACT = 'F').

                     slacpy('Full', N, N, ASAV, LDA, AFAC, LDA );
                     if ( EQUIL || IEQUED > 1 ) {

                        // Compute row and column scale factors to
                        // equilibrate the matrix A.

                        sgeequ(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
                        if ( INFO == 0 && N > 0 ) {
                           if ( lsame( EQUED, 'R' ) ) {
                              ROWCND = ZERO;
                              COLCND = ONE;
                           } else if ( lsame( EQUED, 'C' ) ) {
                              ROWCND = ONE;
                              COLCND = ZERO;
                           } else if ( lsame( EQUED, 'B' ) ) {
                              ROWCND = ZERO;
                              COLCND = ZERO;
                           }

                           // Equilibrate the matrix.

                           slaqge(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                        }
                     }

                     // Save the condition number of the non-equilibrated
                     // system for use in SGET04.

                     if ( EQUIL ) {
                        ROLDO = RCONDO;
                        ROLDI = RCONDI;
                     }

                     // Compute the 1-norm and infinity-norm of A.

                     ANORMO = SLANGE( '1', N, N, AFAC, LDA, RWORK );
                     ANORMI = SLANGE( 'I', N, N, AFAC, LDA, RWORK );

                     // Factor the matrix A.

                    srnamc.SRNAMT = 'SGETRF';
                     sgetrf(N, N, AFAC, LDA, IWORK, INFO );

                     // Form the inverse of A.

                     slacpy('Full', N, N, AFAC, LDA, A, LDA );
                     LWORK = NMAX*max( 3, NRHS );
                    srnamc.SRNAMT = 'SGETRI';
                     sgetri(N, A, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute the 1-norm condition number of A.

                     AINVNM = SLANGE( '1', N, N, A, LDA, RWORK );
                     if ( ANORMO <= ZERO || AINVNM <= ZERO ) {
                        RCONDO = ONE;
                     } else {
                        RCONDO = ( ONE / ANORMO ) / AINVNM;
                     }

                     // Compute the infinity-norm condition number of A.

                     AINVNM = SLANGE( 'I', N, N, A, LDA, RWORK );
                     if ( ANORMI <= ZERO || AINVNM <= ZERO ) {
                        RCONDI = ONE;
                     } else {
                        RCONDI = ( ONE / ANORMI ) / AINVNM;
                     }
                  }

                  for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 50

                     // Do for each value of TRANS.

                     TRANS = TRANSS( ITRAN );
                     if ( ITRAN == 1 ) {
                        RCONDC = RCONDO;
                     } else {
                        RCONDC = RCONDI;
                     }

                     // Restore the matrix A.

                     slacpy('Full', N, N, ASAV, LDA, A, LDA );

                     // Form an exact solution and set the right hand side.

                    srnamc.SRNAMT = 'SLARHS';
                     slarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     slacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT && ITRAN == 1 ) {

                        // --- Test SGESV  ---

                        // Compute the LU factorization of the matrix and
                        // solve the system.

                        slacpy('Full', N, N, A, LDA, AFAC, LDA );
                        slacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'SGESV ';
                        sgesv(N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                        // Check error code from SGESV .

                        if (INFO != IZERO) alaerh( PATH, 'SGESV ', INFO, IZERO, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        // Reconstruct matrix from factors and compute
                        // residual.

                        sget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK, RESULT( 1 ) );
                        NT = 1;
                        if ( IZERO == 0 ) {

                           // Compute residual of the computed solution.

                           slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                           sget02('No transpose', N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                           // Check solution from generated exact solution.

                           sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                           NT = 3;
                        }

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 30
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9999 )'SGESV ', N, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 30
                        NRUN = NRUN + NT;
                     }

                     // --- Test SGESVX ---

                     if ( !PREFAC) slaset( 'Full', N, N, ZERO, ZERO, AFAC, LDA );
                     slaset('Full', N, NRHS, ZERO, ZERO, X, LDA );
                     if ( IEQUED > 1 && N > 0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        slaqge(N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using SGESVX.

                    srnamc.SRNAMT = 'SGESVX';
                     sgesvx(FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                     // Check the error code from SGESVX.

                     if (INFO != IZERO) alaerh( PATH, 'SGESVX', INFO, IZERO, FACT + TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     // Compare WORK(1) from SGESVX with the computed
                     // reciprocal pivot growth factor RPVGRW

                     if ( INFO != 0 && INFO <= N) {
                        RPVGRW = SLANTR( 'M', 'U', 'N', INFO, INFO, AFAC, LDA, WORK );
                        if ( RPVGRW == ZERO ) {
                           RPVGRW = ONE;
                        } else {
                           RPVGRW = SLANGE( 'M', N, INFO, A, LDA, WORK ) / RPVGRW;
                        }
                     } else {
                        RPVGRW = SLANTR( 'M', 'U', 'N', N, N, AFAC, LDA, WORK );
                        if ( RPVGRW == ZERO ) {
                           RPVGRW = ONE;
                        } else {
                           RPVGRW = SLANGE( 'M', N, N, A, LDA, WORK ) / RPVGRW;
                        }
                     }
                     RESULT[7] = ( RPVGRW-WORK( 1 ) ).abs() / max( WORK( 1 ), RPVGRW ) / SLAMCH( 'E' );

                     if ( !PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        sget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1;
                     } else {
                        K1 = 2;
                     }

                     if ( INFO == 0 ) {
                        TRFCON = false;

                        // Compute residual of the computed solution.

                        slacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA );
                        sget02(TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        if( NOFACT || ( PREFAC && lsame( EQUED, 'N' ) ) ) THEN;
                           sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        } else {
                           if ( ITRAN == 1 ) {
                              ROLDC = ROLDO;
                           } else {
                              ROLDC = ROLDI;
                           }
                           sget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        sget07(TRANS, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, true , RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        TRFCON = true;
                     }

                     // Compare RCOND from SGESVX with the computed value
                     // in RCONDC.

                     RESULT[6] = SGET06( RCOND, RCONDC );

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( !TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 40
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K );
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, TRANS, N, IMAT, K, RESULT( K );
                              }
                              NFAIL = NFAIL + 1;
                           }
                        } // 40
                        NRUN = NRUN + NTESTS - K1 + 1;
                     } else {
                        if( RESULT( 1 ) >= THRESH && !PREFAC ) {
                           if( NFAIL == 0 && NERRS == 0 ) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 );
                           } else {
                              WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, TRANS, N, IMAT, 1, RESULT( 1 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }
                        if ( RESULT( 6 ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 );
                           } else {
                              WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, TRANS, N, IMAT, 6, RESULT( 6 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }
                        if ( RESULT( 7 ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'SGESVX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 );
                           } else {
                              WRITE( NOUT, FMT = 9998 )'SGESVX', FACT, TRANS, N, IMAT, 7, RESULT( 7 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }

                     }

                  } // 50
               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${}, N =${.i5}, type ${.i2}, test(${.i2}) =${.g12_5};
 9998 FORMAT(' ${}, FACT=\'${.a1}\', TRANS=\'${.a1}\', N=${.i5}, type ${.i2}, test(${.i1})=${.g12_5};
 9997 FORMAT(' ${}, FACT=\'${.a1}\', TRANS=\'${.a1}\', N=${.i5}, EQUED=\'${.a1}\', type ${.i2}, test(${.i1})=${.g12_5};
      }
