      SUBROUTINE DDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

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
      double             A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), RWORK( * ), S( * ), WORK( * ), X( * ), XACT( * );
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
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, ITRAN, IZERO, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT, N_ERR_BNDS;
      double             AINVNM, AMAX, ANORM, ANORMI, ANORMO, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, RPVGRW_SVXX;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLAMCH, DLANGE, DLANTR, DLA_GERPVGRW;
      // EXTERNAL LSAME, DGET06, DLAMCH, DLANGE, DLANTR, DLA_GERPVGRW
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DERRVX, DGEEQU, DGESV, DGESVX, DGET01, DGET02, DGET04, DGET07, DGETRF, DGETRI, DLACPY, DLAQGE, DLARHS, DLASET, DLATB4, DLATMS, XLAENV, DGESVXX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
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
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               TRANSS / 'N', 'T', 'C' /
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'R', 'C', 'B' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL DERRVX( PATH, NOUT );
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
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 80

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
            if (ZEROT .AND. N.LT.IMAT-4) GO TO 80;

            // Set up parameters with DLATB4 and generate a test matrix
            // with DLATMS.

            dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
            RCONDC = ONE / CNDNUM

            SRNAMT = 'DLATMS'
            dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from DLATMS.

            if ( INFO.NE.0 ) {
               alaerh(PATH, 'DLATMS', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
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
                  dlaset('Full', N, N-IZERO+1, ZERO, ZERO, A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0
            }

            // Save a copy of the matrix A in ASAV.

            dlacpy('Full', N, N, A, LDA, ASAV, LDA );

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
                     if (PREFAC) GO TO 60;
                     RCONDO = ZERO
                     RCONDI = ZERO

                  } else if ( .NOT.NOFACT ) {

                     // Compute the condition number for comparison with
                     // the value returned by DGESVX (FACT = 'N' reuses
                     // the condition number from the previous iteration
                     // with FACT = 'F').

                     dlacpy('Full', N, N, ASAV, LDA, AFAC, LDA );
                     if ( EQUIL .OR. IEQUED.GT.1 ) {

                        // Compute row and column scale factors to
                        // equilibrate the matrix A.

                        dgeequ(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
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

                           dlaqge(N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                        }
                     }

                     // Save the condition number of the non-equilibrated
                     // system for use in DGET04.

                     if ( EQUIL ) {
                        ROLDO = RCONDO
                        ROLDI = RCONDI
                     }

                     // Compute the 1-norm and infinity-norm of A.

                     ANORMO = DLANGE( '1', N, N, AFAC, LDA, RWORK )
                     ANORMI = DLANGE( 'I', N, N, AFAC, LDA, RWORK )

                     // Factor the matrix A.

                     dgetrf(N, N, AFAC, LDA, IWORK, INFO );

                     // Form the inverse of A.

                     dlacpy('Full', N, N, AFAC, LDA, A, LDA );
                     LWORK = NMAX*MAX( 3, NRHS )
                     dgetri(N, A, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute the 1-norm condition number of A.

                     AINVNM = DLANGE( '1', N, N, A, LDA, RWORK )
                     if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDO = ONE
                     } else {
                        RCONDO = ( ONE / ANORMO ) / AINVNM
                     }

                     // Compute the infinity-norm condition number of A.

                     AINVNM = DLANGE( 'I', N, N, A, LDA, RWORK )
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

                     dlacpy('Full', N, N, ASAV, LDA, A, LDA );

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'DLARHS'
                     dlarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C'
                     dlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     if ( NOFACT .AND. ITRAN.EQ.1 ) {

                        // --- Test DGESV  ---

                        // Compute the LU factorization of the matrix and
                        // solve the system.

                        dlacpy('Full', N, N, A, LDA, AFAC, LDA );
                        dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'DGESV '
                        dgesv(N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                        // Check error code from DGESV .

                        if (INFO.NE.IZERO) CALL ALAERH( PATH, 'DGESV ', INFO, IZERO, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        // Reconstruct matrix from factors and compute
                        // residual.

                        dget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK, RESULT( 1 ) );
                        NT = 1
                        if ( IZERO.EQ.0 ) {

                           // Compute residual of the computed solution.

                           dlacpy('Full', N, NRHS, B, LDA, WORK, LDA )                            CALL DGET02( 'No transpose', N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                           // Check solution from generated exact solution.

                           dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                           NT = 3
                        }

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= NT; K++) { // 30
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'DGESV ', N, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1
                           }
                        } // 30
                        NRUN = NRUN + NT
                     }

                     // --- Test DGESVX ---

                     if (.NOT.PREFAC) CALL DLASET( 'Full', N, N, ZERO, ZERO, AFAC, LDA );
                     dlaset('Full', N, NRHS, ZERO, ZERO, X, LDA );
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        dlaqge(N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using DGESVX.

                     SRNAMT = 'DGESVX'
                     dgesvx(FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                     // Check the error code from DGESVX.

                     if (INFO.NE.IZERO) CALL ALAERH( PATH, 'DGESVX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     // Compare WORK(1) from DGESVX with the computed
                     // reciprocal pivot growth factor RPVGRW

                     if ( INFO.NE.0 ) {
                        RPVGRW = DLANTR( 'M', 'U', 'N', INFO, INFO, AFAC, LDA, WORK )
                        if ( RPVGRW.EQ.ZERO ) {
                           RPVGRW = ONE
                        } else {
                           RPVGRW = DLANGE( 'M', N, INFO, A, LDA, WORK ) / RPVGRW
                        }
                     } else {
                        RPVGRW = DLANTR( 'M', 'U', 'N', N, N, AFAC, LDA, WORK )
                        if ( RPVGRW.EQ.ZERO ) {
                           RPVGRW = ONE
                        } else {
                           RPVGRW = DLANGE( 'M', N, N, A, LDA, WORK ) / RPVGRW
                        }
                     }
                     RESULT( 7 ) = ABS( RPVGRW-WORK( 1 ) ) / MAX( WORK( 1 ), RPVGRW ) / DLAMCH( 'E' )

                     if ( .NOT.PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        dget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     if ( INFO.EQ.0 ) {
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        dlacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL DGET02( TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        } else {
                           if ( ITRAN.EQ.1 ) {
                              ROLDC = ROLDO
                           } else {
                              ROLDC = ROLDI
                           }
                           dget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        dget07(TRANS, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, .TRUE., RWORK( NRHS+1 ), RESULT( 4 ) );
                     } else {
                        TRFCON = .TRUE.
                     }

                     // Compare RCOND from DGESVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( .NOT.TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 40
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'DGESVX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'DGESVX', FACT, TRANS, N, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
                        } // 40
                        NRUN = NRUN + 7 - K1
                     } else {
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVX', FACT, TRANS, N, IMAT, 1, RESULT( 1 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 6 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVX', FACT, TRANS, N, IMAT, 6, RESULT( 6 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 7 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVX', FACT, TRANS, N, IMAT, 7, RESULT( 7 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }

                     }

                     // --- Test DGESVXX ---

                     // Restore the matrices A and B.

                     dlacpy('Full', N, N, ASAV, LDA, A, LDA );
                     dlacpy('Full', N, NRHS, BSAV, LDA, B, LDA );
                      if (.NOT.PREFAC) CALL DLASET( 'Full', N, N, ZERO, ZERO, AFAC, LDA );
                     dlaset('Full', N, NRHS, ZERO, ZERO, X, LDA );
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        dlaqge(N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using DGESVXX.

                     SRNAMT = 'DGESVXX'
                     N_ERR_BNDS = 3
                     dgesvxx(FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, IWORK( N+1 ), INFO );

                     // Check the error code from DGESVXX.

                     if (INFO.EQ.N+1) GOTO 50;
                     if ( INFO.NE.IZERO ) {
                        alaerh(PATH, 'DGESVXX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GOTO 50
                     }

                     // Compare rpvgrw_svxx from DGESVXX with the computed
                     // reciprocal pivot growth factor RPVGRW


                     if ( INFO .GT. 0 .AND. INFO .LT. N+1 ) {
                        RPVGRW = DLA_GERPVGRW (N, INFO, A, LDA, AFAC, LDA)
                     } else {
                        RPVGRW = DLA_GERPVGRW (N, N, A, LDA, AFAC, LDA)
                     }
                      RESULT( 7 ) = ABS( RPVGRW-RPVGRW_SVXX ) / MAX( RPVGRW_SVXX, RPVGRW ) / DLAMCH( 'E' )

                     if ( .NOT.PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        dget01(N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     if ( INFO.EQ.0 ) {
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        dlacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL DGET02( TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        } else {
                           if ( ITRAN.EQ.1 ) {
                              ROLDC = ROLDO
                           } else {
                              ROLDC = ROLDI
                           }
                           dget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                        }
                     } else {
                        TRFCON = .TRUE.
                     }

                     // Compare RCOND from DGESVXX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( .NOT.TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 45
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'DGESVXX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'DGESVXX', FACT, TRANS, N, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
                        } // 45
                        NRUN = NRUN + 7 - K1
                     } else {
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVXX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVXX', FACT, TRANS, N, IMAT, 1, RESULT( 1 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 6 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVXX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVXX', FACT, TRANS, N, IMAT, 6, RESULT( 6 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 7 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DGESVXX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DGESVXX', FACT, TRANS, N, IMAT, 7, RESULT( 7 )
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


      // Test Error Bounds from DGESVXX

      debchvxx(THRESH, PATH );

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test(', I2, ') =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', type ', I2, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', EQUED=''', A1, ''', type ', I2, ', test(', I1, ')=', G12.5 )
      RETURN

      // End of DDRVGEX

      }
