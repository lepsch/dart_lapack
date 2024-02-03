      SUBROUTINE ZDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             RWORK( * ), S( * );
      COMPLEX*16         A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, TRFCON, ZEROT;
      String             DIST, EQUED, FACT, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, ITRAN, IZERO, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, AMAX, ANORM, ANORMI, ANORMO, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RDUM( 1 ), RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLAMCH, ZLANGE, ZLANTR;
      // EXTERNAL LSAME, DGET06, DLAMCH, ZLANGE, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGEEQU, ZGESV, ZGESVX, ZGET01, ZGET02, ZGET04, ZGET07, ZGETRF, ZGETRI, ZLACPY, ZLAQGE, ZLARHS, ZLASET, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, MAX
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               TRANSS / 'N', 'T', 'C' /
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'R', 'C', 'B' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 90
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 80

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
            IF( ZEROT .AND. N.LT.IMAT-4 ) GO TO 80

            // Set up parameters with ZLATB4 and generate a test matrix
            // with ZLATMS.

            zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
            RCONDC = ONE / CNDNUM

            SRNAMT = 'ZLATMS'
            zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from ZLATMS.

            if ( INFO.NE.0 ) {
               alaerh(PATH, 'ZLATMS', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 80
            }

            // For types 5-7, zero one or more columns of the matrix to
            // test that INFO is returned correctly.

            if ( ZEROT ) {
               if ( IMAT.EQ.5 ) {
                  IZERO = 1
               } else if ( IMAT.EQ.6 ) {
                  IZERO = N
               } else {
                  IZERO = N / 2 + 1
               }
               IOFF = ( IZERO-1 )*LDA
               if ( IMAT.LT.7 ) {
                  for (I = 1; I <= N; I++) { // 20
                     A( IOFF+I ) = ZERO
                  } // 20
               } else {
                  zlaset('Full', N, N-IZERO+1, DCMPLX( ZERO ), DCMPLX( ZERO ), A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0
            }

            // Save a copy of the matrix A in ASAV.

            zlacpy('Full', N, N, A, LDA, ASAV, LDA );

            for (IEQUED = 1; IEQUED <= 4; IEQUED++) { // 70
               EQUED = EQUEDS( IEQUED )
               if ( IEQUED.EQ.1 ) {
                  NFACT = 3
               } else {
                  NFACT = 1
               }

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 60
                  FACT = FACTS( IFACT )
                  PREFAC = LSAME( FACT, 'F' )
                  NOFACT = LSAME( FACT, 'N' )
                  EQUIL = LSAME( FACT, 'E' )

                  if ( ZEROT ) {
                     IF( PREFAC ) GO TO 60
                     RCONDO = ZERO
                     RCONDI = ZERO

                  } else if ( .NOT.NOFACT ) {

                     // Compute the condition number for comparison with
                     // the value returned by ZGESVX (FACT = 'N' reuses
                     // the condition number from the previous iteration
                     // with FACT = 'F').

                     zlacpy('Full', N, N, ASAV, LDA, AFAC, LDA );
                     if ( EQUIL .OR. IEQUED.GT.1 ) {

                        // Compute row and column scale factors to
                        // equilibrate the matrix A.

                        zgeequ(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
                        if ( INFO.EQ.0 .AND. N.GT.0 ) {
                           if ( LSAME( EQUED, 'R' ) ) {
                              ROWCND = ZERO
                              COLCND = ONE
                           } else if ( LSAME( EQUED, 'C' ) ) {
                              ROWCND = ONE
                              COLCND = ZERO
                           } else if ( LSAME( EQUED, 'B' ) ) {
                              ROWCND = ZERO
                              COLCND = ZERO
                           }

                           // Equilibrate the matrix.

                           zlaqge(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                        }
                     }

                     // Save the condition number of the non-equilibrated
                     // system for use in ZGET04.

                     if ( EQUIL ) {
                        ROLDO = RCONDO
                        ROLDI = RCONDI
                     }

                     // Compute the 1-norm and infinity-norm of A.

                     ANORMO = ZLANGE( '1', N, N, AFAC, LDA, RWORK )
                     ANORMI = ZLANGE( 'I', N, N, AFAC, LDA, RWORK )

                     // Factor the matrix A.

                     SRNAMT = 'ZGETRF'
                     zgetrf(N, N, AFAC, LDA, IWORK, INFO );

                     // Form the inverse of A.

                     zlacpy('Full', N, N, AFAC, LDA, A, LDA );
                     LWORK = NMAX*MAX( 3, NRHS )
                     SRNAMT = 'ZGETRI'
                     zgetri(N, A, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute the 1-norm condition number of A.

                     AINVNM = ZLANGE( '1', N, N, A, LDA, RWORK )
                     if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDO = ONE
                     } else {
                        RCONDO = ( ONE / ANORMO ) / AINVNM
                     }

                     // Compute the infinity-norm condition number of A.

                     AINVNM = ZLANGE( 'I', N, N, A, LDA, RWORK )
                     if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDI = ONE
                     } else {
                        RCONDI = ( ONE / ANORMI ) / AINVNM
                     }
                  }

                  for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 50

                     // Do for each value of TRANS.

                     TRANS = TRANSS( ITRAN )
                     if ( ITRAN.EQ.1 ) {
                        RCONDC = RCONDO
                     } else {
                        RCONDC = RCONDI
                     }

                     // Restore the matrix A.

                     zlacpy('Full', N, N, ASAV, LDA, A, LDA );

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'ZLARHS'
                     zlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C'
                     zlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT .AND. ITRAN.EQ.1 ) {

                        // --- Test ZGESV  ---

                        // Compute the LU factorization of the matrix and
                        // solve the system.

                        zlacpy('Full', N, N, A, LDA, AFAC, LDA );
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'ZGESV '
                        zgesv(N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                        // Check error code from ZGESV .

                        IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'ZGESV ', INFO, IZERO, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        // Reconstruct matrix from factors and compute
                        // residual.

                        zget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK, RESULT( 1 ) );
                        NT = 1
                        if ( IZERO.EQ.0 ) {

                           // Compute residual of the computed solution.

                           zlacpy('Full', N, NRHS, B, LDA, WORK, LDA )                            CALL ZGET02( 'No transpose', N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                           // Check solution from generated exact solution.

                           zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                           NT = 3
                        }

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 30
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'ZGESV ', N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
                        } // 30
                        NRUN = NRUN + NT
                     }

                     // --- Test ZGESVX ---

                     IF( .NOT.PREFAC ) CALL ZLASET( 'Full', N, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, LDA )
                     zlaset('Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA );
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        zlaqge(N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZGESVX.

                     SRNAMT = 'ZGESVX'
                     zgesvx(FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check the error code from ZGESVX.

                     IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'ZGESVX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                     // Compare RWORK(2*NRHS+1) from ZGESVX with the
                     // computed reciprocal pivot growth factor RPVGRW

                     if ( INFO.NE.0 .AND. INFO.LE.N) {
                        RPVGRW = ZLANTR( 'M', 'U', 'N', INFO, INFO, AFAC, LDA, RDUM )
                        if ( RPVGRW.EQ.ZERO ) {
                           RPVGRW = ONE
                        } else {
                           RPVGRW = ZLANGE( 'M', N, INFO, A, LDA, RDUM ) / RPVGRW
                        }
                     } else {
                        RPVGRW = ZLANTR( 'M', 'U', 'N', N, N, AFAC, LDA, RDUM )
                        if ( RPVGRW.EQ.ZERO ) {
                           RPVGRW = ONE
                        } else {
                           RPVGRW = ZLANGE( 'M', N, N, A, LDA, RDUM ) / RPVGRW
                        }
                     }
                     RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) ) / MAX( RWORK( 2*NRHS+1 ), RPVGRW ) / DLAMCH( 'E' )

                     if ( .NOT.PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        zget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     if ( INFO.EQ.0 ) {
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        zlacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL ZGET02( TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        } else {
                           if ( ITRAN.EQ.1 ) {
                              ROLDC = ROLDO
                           } else {
                              ROLDC = ROLDI
                           }
                           zget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        zget07(TRANS, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, .TRUE., RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        TRFCON = .TRUE.
                     }

                     // Compare RCOND from ZGESVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( .NOT.TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 40
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'ZGESVX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'ZGESVX', FACT, TRANS, N, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
                        } // 40
                        NRUN = NRUN + NTESTS - K1 + 1
                     } else {
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'ZGESVX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'ZGESVX', FACT, TRANS, N, IMAT, 1, RESULT( 1 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 6 ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'ZGESVX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'ZGESVX', FACT, TRANS, N, IMAT, 6, RESULT( 6 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 7 ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'ZGESVX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'ZGESVX', FACT, TRANS, N, IMAT, 7, RESULT( 7 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }

                     }

                  } // 50
               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test(', I2, ') =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', type ', I2, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', EQUED=''', A1, ''', type ', I2, ', test(', I1, ')=', G12.5 )
      RETURN

      // End of ZDRVGE

      }
