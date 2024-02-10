      void zdrvgb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA, AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, final Array<int> IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                LA, LAFB, NN, NOUT, NRHS;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             RWORK( * ), S( * );
      Complex         A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      bool               EQUIL, NOFACT, PREFAC, TRFCON, ZEROT;
      String             DIST, EQUED, FACT, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, I1, I2, IEQUED, IFACT, IKL, IKU, IMAT, IN, INFO, IOFF, ITRAN, IZERO, J, K, K1, KL, KU, LDA, LDAFB, LDB, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NKL, NKU, NRUN, NT, N_ERR_BNDS;
      double             AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, RPVGRW_SVXX;
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DGET06, DLAMCH, ZLANGB, ZLANGE, ZLANTB, ZLA_GBRPVGRW       EXTERNAL           lsame, DGET06, DLAMCH, ZLANGB, ZLANGE, ZLANTB, ZLA_GBRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGBEQU, ZGBSV, ZGBSVX, ZGBT01, ZGBT02, ZGBT05, ZGBTRF, ZGBTRS, ZGET04, ZLACPY, ZLAQGB, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZGBSVXX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, MAX, MIN
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
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const TRANSS = [ 'N', 'T', 'C' ];
      const FACTS = [ 'F', 'N', 'E' ];
      const EQUEDS = [ 'N', 'R', 'C', 'B' ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'GB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrvx( PATH, NOUT );
      infoc.INFOT = 0;

      // Set the block size and minimum block size for testing.

      NB = 1;
      NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 150
         N = NVAL( IN );
         LDB = max( N, 1 );
         XTYPE = 'N';

         // Set limits on the number of loop iterations.

         NKL = max( 1, min( N, 4 ) );
         if (N == 0) NKL = 1;
         NKU = NKL;
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IKL = 1; IKL <= NKL; IKL++) { // 140

            // Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
            // it easier to skip redundant values for small values of N.

            if ( IKL == 1 ) {
               KL = 0;
            } else if ( IKL == 2 ) {
               KL = max( N-1, 0 );
            } else if ( IKL == 3 ) {
               KL = ( 3*N-1 ) / 4;
            } else if ( IKL == 4 ) {
               KL = ( N+1 ) / 4;
            }
            for (IKU = 1; IKU <= NKU; IKU++) { // 130

               // Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
               // makes it easier to skip redundant values for small
               // values of N.

               if ( IKU == 1 ) {
                  KU = 0;
               } else if ( IKU == 2 ) {
                  KU = max( N-1, 0 );
               } else if ( IKU == 3 ) {
                  KU = ( 3*N-1 ) / 4;
               } else if ( IKU == 4 ) {
                  KU = ( N+1 ) / 4;
               }

               // Check that A and AFB are big enough to generate this
               // matrix.

               LDA = KL + KU + 1;
               LDAFB = 2*KL + KU + 1;
               if ( LDA*N > LA || LDAFB*N > LAFB ) {
                  if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                  if ( LDA*N > LA ) {
                     WRITE( NOUT, FMT = 9999 )LA, N, KL, KU, N*( KL+KU+1 );
                     NERRS = NERRS + 1;
                  }
                  if ( LDAFB*N > LAFB ) {
                     WRITE( NOUT, FMT = 9998 )LAFB, N, KL, KU, N*( 2*KL+KU+1 );
                     NERRS = NERRS + 1;
                  }
                  GO TO 130;
               }

               for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 120

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  if( !DOTYPE( IMAT ) ) GO TO 120;

                  // Skip types 2, 3, or 4 if the matrix is too small.

                  ZEROT = IMAT >= 2 && IMAT <= 4;
                  if (ZEROT && N < IMAT-1) GO TO 120;

                  // Set up parameters with ZLATB4 and generate a
                  // test matrix with ZLATMS.

                  zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
                  RCONDC = ONE / CNDNUM;

                 srnamc.SRNAMT = 'ZLATMS';
                  zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK, INFO );

                  // Check the error code from ZLATMS.

                  if ( INFO != 0 ) {
                     alaerh(PATH, 'ZLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 120;
                  }

                  // For types 2, 3, and 4, zero one or more columns of
                  // the matrix to test that INFO is returned correctly.

                  IZERO = 0;
                  if ( ZEROT ) {
                     if ( IMAT == 2 ) {
                        IZERO = 1;
                     } else if ( IMAT == 3 ) {
                        IZERO = N;
                     } else {
                        IZERO = N / 2 + 1;
                     }
                     IOFF = ( IZERO-1 )*LDA;
                     if ( IMAT < 4 ) {
                        I1 = max( 1, KU+2-IZERO );
                        I2 = min( KL+KU+1, KU+1+( N-IZERO ) );
                        for (I = I1; I <= I2; I++) { // 20
                           A[IOFF+I] = ZERO;
                        } // 20
                     } else {
                        for (J = IZERO; J <= N; J++) { // 40
                           DO 30 I = max( 1, KU+2-J ), min( KL+KU+1, KU+1+( N-J ) );
                              A[IOFF+I] = ZERO;
                           } // 30
                           IOFF = IOFF + LDA;
                        } // 40
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  zlacpy('Full', KL+KU+1, N, A, LDA, ASAV, LDA );

                  for (IEQUED = 1; IEQUED <= 4; IEQUED++) { // 110
                     EQUED = EQUEDS( IEQUED );
                     if ( IEQUED == 1 ) {
                        NFACT = 3;
                     } else {
                        NFACT = 1;
                     }

                     for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 100
                        FACT = FACTS( IFACT );
                        PREFAC = lsame( FACT, 'F' );
                        NOFACT = lsame( FACT, 'N' );
                        EQUIL = lsame( FACT, 'E' );

                        if ( ZEROT ) {
                           if (PREFAC) GO TO 100;
                           RCONDO = ZERO;
                           RCONDI = ZERO;

                        } else if ( !NOFACT ) {

                           // Compute the condition number for comparison
                           // with the value returned by DGESVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           zlacpy('Full', KL+KU+1, N, ASAV, LDA, AFB( KL+1 ), LDAFB );
                           if ( EQUIL || IEQUED > 1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              zgbequ(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
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

                                 zlaqgb(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in ZGET04.

                           if ( EQUIL ) {
                              ROLDO = RCONDO;
                              ROLDI = RCONDI;
                           }

                           // Compute the 1-norm and infinity-norm of A.

                           ANORMO = ZLANGB( '1', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )                            ANORMI = ZLANGB( 'I', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK );

                           // Factor the matrix A.

                           zgbtrf(N, N, KL, KU, AFB, LDAFB, IWORK, INFO );

                           // Form the inverse of A.

                           zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, LDB );
                          srnamc.SRNAMT = 'ZGBTRS';
                           zgbtrs('No transpose', N, KL, KU, N, AFB, LDAFB, IWORK, WORK, LDB, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = ZLANGE( '1', N, N, WORK, LDB, RWORK );
                           if ( ANORMO <= ZERO || AINVNM <= ZERO ) {
                              RCONDO = ONE;
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM;
                           }

                           // Compute the infinity-norm condition number
                           // of A.

                           AINVNM = ZLANGE( 'I', N, N, WORK, LDB, RWORK );
                           if ( ANORMI <= ZERO || AINVNM <= ZERO ) {
                              RCONDI = ONE;
                           } else {
                              RCONDI = ( ONE / ANORMI ) / AINVNM;
                           }
                        }

                        for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 90

                           // Do for each value of TRANS.

                           TRANS = TRANSS( ITRAN );
                           if ( ITRAN == 1 ) {
                              RCONDC = RCONDO;
                           } else {
                              RCONDC = RCONDI;
                           }

                           // Restore the matrix A.

                           zlacpy('Full', KL+KU+1, N, ASAV, LDA, A, LDA );

                           // Form an exact solution and set the right hand
                           // side.

                          srnamc.SRNAMT = 'ZLARHS';
                           zlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO );
                           XTYPE = 'C';
                           zlacpy('Full', N, NRHS, B, LDB, BSAV, LDB );

                           if ( NOFACT && ITRAN == 1 ) {

                              // --- Test ZGBSV  ---

                              // Compute the LU factorization of the matrix
                              // and solve the system.

                              zlacpy('Full', KL+KU+1, N, A, LDA, AFB( KL+1 ), LDAFB );
                              zlacpy('Full', N, NRHS, B, LDB, X, LDB );

                             srnamc.SRNAMT = 'ZGBSV ';
                              zgbsv(N, KL, KU, NRHS, AFB, LDAFB, IWORK, X, LDB, INFO );

                              // Check error code from ZGBSV .

                              if (INFO != IZERO) alaerh( PATH, 'ZGBSV ', INFO, IZERO, ' ', N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                              // Reconstruct matrix from factors and
                              // compute residual.

                              zgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              NT = 1;
                              if ( IZERO == 0 ) {

                                 // Compute residual of the computed
                                 // solution.

                                 zlacpy('Full', N, NRHS, B, LDB, WORK, LDB );
                                 zgbt02('No transpose', N, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

                                 // Check solution from generated exact
                                 // solution.

                                 zget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                                 NT = 3;
                              }

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (K = 1; K <= NT; K++) { // 50
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9997 )'ZGBSV ', N, KL, KU, IMAT, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              } // 50
                              NRUN = NRUN + NT;
                           }

                           // --- Test ZGBSVX ---

                           if ( !PREFAC) zlaset( 'Full', 2*KL+KU+1, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFB, LDAFB );
                           zlaset('Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDB );
                           if ( IEQUED > 1 && N > 0 ) {

                              // Equilibrate the matrix if FACT = 'F' and
                              // EQUED = 'R', 'C', or 'B'.

                              zlaqgb(N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                           }

                           // Solve the system and compute the condition
                           // number and error bounds using ZGBSVX.

                          srnamc.SRNAMT = 'ZGBSVX';
                           zgbsvx(FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( LDB+1 ), B, LDB, X, LDB, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                           // Check the error code from ZGBSVX.

                           if (INFO != IZERO) alaerh( PATH, 'ZGBSVX', INFO, IZERO, FACT + TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                           // Compare RWORK(2*NRHS+1) from ZGBSVX with the
                           // computed reciprocal pivot growth RPVGRW

                           if ( INFO != 0 ) {
                              ANRMPV = ZERO;
                              for (J = 1; J <= INFO; J++) { // 70
                                 for (I = max( KU+2-J, 1 ); I <= min( N+KU+1-J, KL+KU+1 ); I++) { // 60
                                    ANRMPV = max( ANRMPV, ABS( A( I+( J-1 )*LDA ) ) )
                                 } // 60
                              } // 70
                              RPVGRW = ZLANTB( 'M', 'U', 'N', INFO, min( INFO-1, KL+KU ), AFB( max( 1, KL+KU+2-INFO ) ), LDAFB, RDUM );
                              if ( RPVGRW == ZERO ) {
                                 RPVGRW = ONE;
                              } else {
                                 RPVGRW = ANRMPV / RPVGRW;
                              }
                           } else {
                              RPVGRW = ZLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, RDUM );
                              if ( RPVGRW == ZERO ) {
                                 RPVGRW = ONE;
                              } else {
                                 RPVGRW = ZLANGB( 'M', N, KL, KU, A, LDA, RDUM ) / RPVGRW;
                              }
                           }
                           RESULT[7] = ( RPVGRW-RWORK( 2*NRHS+1 ) ).abs() / max( RWORK( 2*NRHS+1 ), RPVGRW ) / dlamch( 'E' );

                           if ( !PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              zgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              K1 = 1;
                           } else {
                              K1 = 2;
                           }

                           if ( INFO == 0 ) {
                              TRFCON = false;

                              // Compute residual of the computed solution.

                              zlacpy('Full', N, NRHS, BSAV, LDB, WORK, LDB );
                              zgbt02(TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                              // Check solution from generated exact
                              // solution.

                              if( NOFACT || ( PREFAC && lsame( EQUED, 'N' ) ) ) THEN;
                                 zget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                              } else {
                                 if ( ITRAN == 1 ) {
                                    ROLDC = ROLDO;
                                 } else {
                                    ROLDC = ROLDI;
                                 }
                                 zget04(N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) );
                              }

                              // Check the error bounds from iterative
                              // refinement.

                              zgbt05(TRANS, N, KL, KU, NRHS, ASAV, LDA, BSAV, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                           } else {
                              TRFCON = true;
                           }

                           // Compare RCOND from ZGBSVX with the computed
                           // value in RCONDC.

                           RESULT[6] = DGET06( RCOND, RCONDC );

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( !TRFCON ) {
                              for (K = K1; K <= NTESTS; K++) { // 80
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                                    if ( PREFAC ) {
                                       WRITE( NOUT, FMT = 9995 ) 'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K );
                                    } else {
                                       WRITE( NOUT, FMT = 9996 ) 'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K );
                                    }
                                    NFAIL = NFAIL + 1;
                                 }
                              } // 80
                              NRUN = NRUN + 7 - K1;
                           } else {
                              if( RESULT( 1 ) >= THRESH && !PREFAC ) {
                                 if( NFAIL == 0 && NERRS == 0 ) aladhd( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 );
                                 }
                                 NFAIL = NFAIL + 1;
                                 NRUN = NRUN + 1;
                              }
                              if ( RESULT( 6 ) >= THRESH ) {
                                 if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 );
                                 }
                                 NFAIL = NFAIL + 1;
                                 NRUN = NRUN + 1;
                              }
                              if ( RESULT( 7 ) >= THRESH ) {
                                 if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 );
                                 }
                                 NFAIL = NFAIL + 1;
                                 NRUN = NRUN + 1;
                              }
                           }

                     // --- Test ZGBSVXX ---

                     // Restore the matrices A and B.

                      // write(*,*) 'begin zgbsvxx testing'
                      zlacpy('Full', KL+KU+1, N, ASAV, LDA, A, LDA );
                     zlacpy('Full', N, NRHS, BSAV, LDB, B, LDB );
                      if ( !PREFAC) zlaset( 'Full', 2*KL+KU+1, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFB, LDAFB );
                     zlaset('Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDB );
                     if ( IEQUED > 1 && N > 0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        zlaqgb(N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZGBSVXX.

                    srnamc.SRNAMT = 'ZGBSVXX';
                     N_ERR_BNDS = 3;
                     zgbsvxx(FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB, X, LDB, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, RWORK, INFO );

                     // Check the error code from ZGBSVXX.

                     if (INFO == N+1) GOTO 90;
                     if ( INFO != IZERO ) {
                        alaerh(PATH, 'ZGBSVXX', INFO, IZERO, FACT + TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GOTO 90;
                     }

                     // Compare rpvgrw_svxx from ZGESVXX with the computed
                     // reciprocal pivot growth factor RPVGRW


                     if ( INFO > 0 && INFO < N+1 ) {
                        RPVGRW = ZLA_GBRPVGRW(N, KL, KU, INFO, A, LDA, AFB, LDAFB);
                     } else {
                        RPVGRW = ZLA_GBRPVGRW(N, KL, KU, N, A, LDA, AFB, LDAFB);
                     }
                      RESULT[7] = ( RPVGRW-rpvgrw_svxx ).abs() / max( rpvgrw_svxx, RPVGRW ) / dlamch( 'E' );

                     if ( !PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        zgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1;
                     } else {
                        K1 = 2;
                     }

                     if ( INFO == 0 ) {
                        TRFCON = false;

                        // Compute residual of the computed solution.

                        zlacpy('Full', N, NRHS, BSAV, LDB, WORK, LDB );
                        zgbt02(TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        if( NOFACT || ( PREFAC && lsame( EQUED, 'N' ) ) ) THEN;
                           zget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                        } else {
                           if ( ITRAN == 1 ) {
                              ROLDC = ROLDO;
                           } else {
                              ROLDC = ROLDI;
                           }
                           zget04(N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) );
                        }
                     } else {
                        TRFCON = true;
                     }

                     // Compare RCOND from ZGBSVXX with the computed value
                     // in RCONDC.

                     RESULT[6] = DGET06( RCOND, RCONDC );

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( !TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 45
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K );
                              } else {
                                 WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K );
                              }
                              NFAIL = NFAIL + 1;
                           }
                        } // 45
                        NRUN = NRUN + 7 - K1;
                     } else {
                        if( RESULT( 1 ) >= THRESH && !PREFAC ) {
                           if( NFAIL == 0 && NERRS == 0 ) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 );
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }
                        if ( RESULT( 6 ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 );
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }
                        if ( RESULT( 7 ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 );
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 );
                           }
                           NFAIL = NFAIL + 1;
                           NRUN = NRUN + 1;
                        }

                     }

                        } // 90
                     } // 100
                  } // 110
               } // 120
            } // 130
         } // 140
      } // 150

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );


      // Test Error Bounds from ZGBSVXX

      zebchvxx(THRESH, PATH);

 9999 FORMAT( ' *** In ZDRVGB, LA=${.i5} is too small for N=${.i5}, KU=${.i5}, KL=', I5, / ' ==> Increase LA to at least ${.i5}');
 9998 FORMAT( ' *** In ZDRVGB, LAFB=${.i5} is too small for N=${.i5}, KU=${.i5}, KL=', I5, / ' ==> Increase LAFB to at least ${.i5}');
 9997 FORMAT(' ${}, N=${.i5}, KL=${.i5}, KU=${.i5}, type ${.i1}, test(${.i1})=${.g12_5};
 9996 FORMAT(' ${}( ''${.a1}'',''${.a1}'',${.i5},${.i5},${.i5},...), type ${.i1}, test(${.i1})=${.g12_5};
 9995 FORMAT(' ${}( ''${.a1}'',''${.a1}'',${.i5},${.i5},${.i5},...), EQUED=''${.a1}'', type ${.i1}, test(${.i1})=${.g12_5};

      }
