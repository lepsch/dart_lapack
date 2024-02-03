      SUBROUTINE DDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA, AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                LA, LAFB, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ), RWORK( * ), S( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, TRFCON, ZEROT;
      String             DIST, EQUED, FACT, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, I1, I2, IEQUED, IFACT, IKL, IKU, IMAT, IN, INFO, IOFF, ITRAN, IZERO, J, K, K1, KL, KU, LDA, LDAFB, LDB, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NKL, NKU, NRUN, NT;
      double             AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLAMCH, DLANGB, DLANGE, DLANTB;
      // EXTERNAL LSAME, DGET06, DLAMCH, DLANGB, DLANGE, DLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DERRVX, DGBEQU, DGBSV, DGBSVX, DGBT01, DGBT02, DGBT05, DGBTRF, DGBTRS, DGET04, DLACPY, DLAQGB, DLARHS, DLASET, DLATB4, DLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
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
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /;
      DATA               TRANSS / 'N', 'T', 'C' /;
      DATA               FACTS / 'F', 'N', 'E' /;
      DATA               EQUEDS / 'N', 'R', 'C', 'B' /;
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'GB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) CALL DERRVX( PATH, NOUT );
      INFOT = 0;

      // Set the block size and minimum block size for testing.

      NB = 1;
      NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 150
         N = NVAL( IN );
         LDB = MAX( N, 1 );
         XTYPE = 'N';

         // Set limits on the number of loop iterations.

         NKL = MAX( 1, MIN( N, 4 ) );
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
               KL = MAX( N-1, 0 );
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
                  KU = MAX( N-1, 0 );
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
                  if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
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

                  IF( !DOTYPE( IMAT ) ) GO TO 120;

                  // Skip types 2, 3, or 4 if the matrix is too small.

                  ZEROT = IMAT >= 2 && IMAT <= 4;
                  if (ZEROT && N < IMAT-1) GO TO 120;

                  // Set up parameters with DLATB4 and generate a
                  // test matrix with DLATMS.

                  dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
                  RCONDC = ONE / CNDNUM;

                  SRNAMT = 'DLATMS';
                  dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK, INFO );

                  // Check the error code from DLATMS.

                  if ( INFO != 0 ) {
                     alaerh(PATH, 'DLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
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
                        I1 = MAX( 1, KU+2-IZERO );
                        I2 = MIN( KL+KU+1, KU+1+( N-IZERO ) );
                        for (I = I1; I <= I2; I++) { // 20
                           A( IOFF+I ) = ZERO;
                        } // 20
                     } else {
                        for (J = IZERO; J <= N; J++) { // 40
                           DO 30 I = MAX( 1, KU+2-J ), MIN( KL+KU+1, KU+1+( N-J ) );
                              A( IOFF+I ) = ZERO;
                           } // 30
                           IOFF = IOFF + LDA;
                        } // 40
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  dlacpy('Full', KL+KU+1, N, A, LDA, ASAV, LDA );

                  for (IEQUED = 1; IEQUED <= 4; IEQUED++) { // 110
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
                           RCONDO = ZERO;
                           RCONDI = ZERO;

                        } else if ( !NOFACT ) {

                           // Compute the condition number for comparison
                           // with the value returned by DGESVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           dlacpy('Full', KL+KU+1, N, ASAV, LDA, AFB( KL+1 ), LDAFB );
                           if ( EQUIL || IEQUED > 1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              dgbequ(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
                              if ( INFO == 0 && N > 0 ) {
                                 if ( LSAME( EQUED, 'R' ) ) {
                                    ROWCND = ZERO;
                                    COLCND = ONE;
                                 } else if ( LSAME( EQUED, 'C' ) ) {
                                    ROWCND = ONE;
                                    COLCND = ZERO;
                                 } else if ( LSAME( EQUED, 'B' ) ) {
                                    ROWCND = ZERO;
                                    COLCND = ZERO;
                                 }

                                 // Equilibrate the matrix.

                                 dlaqgb(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in DGET04.

                           if ( EQUIL ) {
                              ROLDO = RCONDO;
                              ROLDI = RCONDI;
                           }

                           // Compute the 1-norm and infinity-norm of A.

                           ANORMO = DLANGB( '1', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )                            ANORMI = DLANGB( 'I', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK );

                           // Factor the matrix A.

                           dgbtrf(N, N, KL, KU, AFB, LDAFB, IWORK, INFO );

                           // Form the inverse of A.

                           dlaset('Full', N, N, ZERO, ONE, WORK, LDB );
                           SRNAMT = 'DGBTRS';
                           dgbtrs('No transpose', N, KL, KU, N, AFB, LDAFB, IWORK, WORK, LDB, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = DLANGE( '1', N, N, WORK, LDB, RWORK );
                           if ( ANORMO <= ZERO || AINVNM <= ZERO ) {
                              RCONDO = ONE;
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM;
                           }

                           // Compute the infinity-norm condition number
                           // of A.

                           AINVNM = DLANGE( 'I', N, N, WORK, LDB, RWORK );
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

                           dlacpy('Full', KL+KU+1, N, ASAV, LDA, A, LDA );

                           // Form an exact solution and set the right hand
                           // side.

                           SRNAMT = 'DLARHS';
                           dlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO );
                           XTYPE = 'C';
                           dlacpy('Full', N, NRHS, B, LDB, BSAV, LDB );

                           if ( NOFACT && ITRAN == 1 ) {

                              // --- Test DGBSV  ---

                              // Compute the LU factorization of the matrix
                              // and solve the system.

                              dlacpy('Full', KL+KU+1, N, A, LDA, AFB( KL+1 ), LDAFB );
                              dlacpy('Full', N, NRHS, B, LDB, X, LDB );

                              SRNAMT = 'DGBSV ';
                              dgbsv(N, KL, KU, NRHS, AFB, LDAFB, IWORK, X, LDB, INFO );

                              // Check error code from DGBSV .

                              if (INFO != IZERO) CALL ALAERH( PATH, 'DGBSV ', INFO, IZERO, ' ', N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                              // Reconstruct matrix from factors and
                              // compute residual.

                              dgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              NT = 1;
                              if ( IZERO == 0 ) {

                                 // Compute residual of the computed
                                 // solution.

                                 dlacpy('Full', N, NRHS, B, LDB, WORK, LDB );
                                 dgbt02('No transpose', N, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

                                 // Check solution from generated exact
                                 // solution.

                                 dget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                                 NT = 3;
                              }

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (K = 1; K <= NT; K++) { // 50
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9997 )'DGBSV ', N, KL, KU, IMAT, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              } // 50
                              NRUN = NRUN + NT;
                           }

                           // --- Test DGBSVX ---

                           if ( !PREFAC) CALL DLASET( 'Full', 2*KL+KU+1, N, ZERO, ZERO, AFB, LDAFB );
                           dlaset('Full', N, NRHS, ZERO, ZERO, X, LDB );
                           if ( IEQUED > 1 && N > 0 ) {

                              // Equilibrate the matrix if FACT = 'F' and
                              // EQUED = 'R', 'C', or 'B'.

                              dlaqgb(N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                           }

                           // Solve the system and compute the condition
                           // number and error bounds using DGBSVX.

                           SRNAMT = 'DGBSVX';
                           dgbsvx(FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB, X, LDB, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                           // Check the error code from DGBSVX.

                           if (INFO != IZERO) CALL ALAERH( PATH, 'DGBSVX', INFO, IZERO, FACT // TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                           // Compare WORK(1) from DGBSVX with the computed
                           // reciprocal pivot growth factor RPVGRW

                           if ( INFO != 0 && INFO <= N) {
                              ANRMPV = ZERO;
                              for (J = 1; J <= INFO; J++) { // 70
                                 DO 60 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )                                     ANRMPV = MAX( ANRMPV, ABS( A( I+( J-1 )*LDA ) ) );
                                 } // 60
                              } // 70
                              RPVGRW = DLANTB( 'M', 'U', 'N', INFO, MIN( INFO-1, KL+KU ), AFB( MAX( 1, KL+KU+2-INFO ) ), LDAFB, WORK );
                              if ( RPVGRW == ZERO ) {
                                 RPVGRW = ONE;
                              } else {
                                 RPVGRW = ANRMPV / RPVGRW;
                              }
                           } else {
                              RPVGRW = DLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, WORK );
                              if ( RPVGRW == ZERO ) {
                                 RPVGRW = ONE;
                              } else {
                                 RPVGRW = DLANGB( 'M', N, KL, KU, A, LDA, WORK ) / RPVGRW;
                              }
                           }
                           RESULT( 7 ) = ABS( RPVGRW-WORK( 1 ) ) / MAX( WORK( 1 ), RPVGRW ) / DLAMCH( 'E' );

                           if ( !PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              dgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              K1 = 1;
                           } else {
                              K1 = 2;
                           }

                           if ( INFO == 0 ) {
                              TRFCON = false;

                              // Compute residual of the computed solution.

                              dlacpy('Full', N, NRHS, BSAV, LDB, WORK, LDB );
                              dgbt02(TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                              // Check solution from generated exact
                              // solution.

                              IF( NOFACT || ( PREFAC && LSAME( EQUED, 'N' ) ) ) THEN;
                                 dget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                              } else {
                                 if ( ITRAN == 1 ) {
                                    ROLDC = ROLDO;
                                 } else {
                                    ROLDC = ROLDI;
                                 }
                                 dget04(N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) );
                              }

                              // Check the error bounds from iterative
                              // refinement.

                              dgbt05(TRANS, N, KL, KU, NRHS, ASAV, LDA, B, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                           } else {
                              TRFCON = true;
                           }

                           // Compare RCOND from DGBSVX with the computed
                           // value in RCONDC.

                           RESULT( 6 ) = DGET06( RCOND, RCONDC );

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( !TRFCON ) {
                              for (K = K1; K <= NTESTS; K++) { // 80
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                                    if ( PREFAC ) {
                                       WRITE( NOUT, FMT = 9995 ) 'DGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K );
                                    } else {
                                       WRITE( NOUT, FMT = 9996 ) 'DGBSVX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K );
                                    }
                                    NFAIL = NFAIL + 1;
                                 }
                              } // 80
                              NRUN = NRUN + NTESTS - K1 + 1;
                           } else {
                              IF( RESULT( 1 ) >= THRESH && !PREFAC ) THEN                                  IF( NFAIL == 0 && NERRS == 0 ) CALL ALADHD( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'DGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'DGBSVX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 );
                                 }
                                 NFAIL = NFAIL + 1;
                                 NRUN = NRUN + 1;
                              }
                              if ( RESULT( 6 ) >= THRESH ) {
                                 if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'DGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'DGBSVX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 );
                                 }
                                 NFAIL = NFAIL + 1;
                                 NRUN = NRUN + 1;
                              }
                              if ( RESULT( 7 ) >= THRESH ) {
                                 if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'DGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 );
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'DGBSVX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 );
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

 9999 FORMAT( ' *** In DDRVGB, LA=', I5, ' is too small for N=', I5, ', KU=', I5, ', KL=', I5, / ' ==> Increase LA to at least ', I5 );
 9998 FORMAT( ' *** In DDRVGB, LAFB=', I5, ' is too small for N=', I5, ', KU=', I5, ', KL=', I5, / ' ==> Increase LAFB to at least ', I5 );
 9997 FORMAT( 1X, A, ', N=', I5, ', KL=', I5, ', KU=', I5, ', type ', I1, ', test(', I1, ')=', G12.5 );
 9996 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', I5, ',...), type ', I1, ', test(', I1, ')=', G12.5 );
 9995 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', I5, ',...), EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=', G12.5 );

      return;

      // End of DDRVGB

      }
