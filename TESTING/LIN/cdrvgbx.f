      SUBROUTINE CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA, AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                LA, LAFB, NN, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
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
      int                I, I1, I2, IEQUED, IFACT, IKL, IKU, IMAT, IN, INFO, IOFF, ITRAN, IZERO, J, K, K1, KL, KU, LDA, LDAFB, LDB, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NKL, NKU, NRUN, NT, N_ERR_BNDS;
      REAL               AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, RPVGRW_SVXX
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS,3 ), ERRBNDS_C( NRHS, 3 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGB, CLANGE, CLANTB, SGET06, SLAMCH, CLA_GBRPVGRW       EXTERNAL           LSAME, CLANGB, CLANGE, CLANTB, SGET06, SLAMCH, CLA_GBRPVGRW
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CERRVX, CGBEQU, CGBSV, CGBSVX, CGBT01, CGBT02, CGBT05, CGBTRF, CGBTRS, CGET04, CLACPY, CLAQGB, CLARHS, CLASET, CLATB4, CLATMS, XLAENV, CGBSVXX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, MIN
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

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'GB'
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

      for (IN = 1; IN <= NN; IN++) { // 150
         N = NVAL( IN )
         LDB = MAX( N, 1 )
         XTYPE = 'N'

         // Set limits on the number of loop iterations.

         NKL = MAX( 1, MIN( N, 4 ) )
         if (N.EQ.0) NKL = 1;
         NKU = NKL
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IKL = 1; IKL <= NKL; IKL++) { // 140

            // Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
            // it easier to skip redundant values for small values of N.

            if ( IKL.EQ.1 ) {
               KL = 0
            } else if ( IKL.EQ.2 ) {
               KL = MAX( N-1, 0 )
            } else if ( IKL.EQ.3 ) {
               KL = ( 3*N-1 ) / 4
            } else if ( IKL.EQ.4 ) {
               KL = ( N+1 ) / 4
            }
            for (IKU = 1; IKU <= NKU; IKU++) { // 130

               // Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
               // makes it easier to skip redundant values for small
               // values of N.

               if ( IKU.EQ.1 ) {
                  KU = 0
               } else if ( IKU.EQ.2 ) {
                  KU = MAX( N-1, 0 )
               } else if ( IKU.EQ.3 ) {
                  KU = ( 3*N-1 ) / 4
               } else if ( IKU.EQ.4 ) {
                  KU = ( N+1 ) / 4
               }

               // Check that A and AFB are big enough to generate this
               // matrix.

               LDA = KL + KU + 1
               LDAFB = 2*KL + KU + 1
               if ( LDA*N.GT.LA .OR. LDAFB*N.GT.LAFB ) {
                  if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                  if ( LDA*N.GT.LA ) {
                     WRITE( NOUT, FMT = 9999 )LA, N, KL, KU, N*( KL+KU+1 )
                     NERRS = NERRS + 1
                  }
                  if ( LDAFB*N.GT.LAFB ) {
                     WRITE( NOUT, FMT = 9998 )LAFB, N, KL, KU, N*( 2*KL+KU+1 )
                     NERRS = NERRS + 1
                  }
                  GO TO 130
               }

               for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 120

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

                  // Skip types 2, 3, or 4 if the matrix is too small.

                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  if (ZEROT .AND. N.LT.IMAT-1) GO TO 120;

                  // Set up parameters with CLATB4 and generate a
                  // test matrix with CLATMS.

                  clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
                  RCONDC = ONE / CNDNUM

                  SRNAMT = 'CLATMS'
                  clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK, INFO );

                  // Check the error code from CLATMS.

                  if ( INFO.NE.0 ) {
                     alaerh(PATH, 'CLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 120
                  }

                  // For types 2, 3, and 4, zero one or more columns of
                  // the matrix to test that INFO is returned correctly.

                  IZERO = 0
                  if ( ZEROT ) {
                     if ( IMAT.EQ.2 ) {
                        IZERO = 1
                     } else if ( IMAT.EQ.3 ) {
                        IZERO = N
                     } else {
                        IZERO = N / 2 + 1
                     }
                     IOFF = ( IZERO-1 )*LDA
                     if ( IMAT.LT.4 ) {
                        I1 = MAX( 1, KU+2-IZERO )
                        I2 = MIN( KL+KU+1, KU+1+( N-IZERO ) )
                        for (I = I1; I <= I2; I++) { // 20
                           A( IOFF+I ) = ZERO
                        } // 20
                     } else {
                        for (J = IZERO; J <= N; J++) { // 40
                           DO 30 I = MAX( 1, KU+2-J ), MIN( KL+KU+1, KU+1+( N-J ) )
                              A( IOFF+I ) = ZERO
                           } // 30
                           IOFF = IOFF + LDA
                        } // 40
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  clacpy('Full', KL+KU+1, N, A, LDA, ASAV, LDA );

                  for (IEQUED = 1; IEQUED <= 4; IEQUED++) { // 110
                     EQUED = EQUEDS( IEQUED )
                     if ( IEQUED.EQ.1 ) {
                        NFACT = 3
                     } else {
                        NFACT = 1
                     }

                     for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 100
                        FACT = FACTS( IFACT )
                        PREFAC = LSAME( FACT, 'F' )
                        NOFACT = LSAME( FACT, 'N' )
                        EQUIL = LSAME( FACT, 'E' )

                        if ( ZEROT ) {
                           if (PREFAC) GO TO 100;
                           RCONDO = ZERO
                           RCONDI = ZERO

                        } else if ( .NOT.NOFACT ) {

                           // Compute the condition number for comparison
                           // with the value returned by SGESVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           clacpy('Full', KL+KU+1, N, ASAV, LDA, AFB( KL+1 ), LDAFB );
                           if ( EQUIL .OR. IEQUED.GT.1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              cgbequ(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO );
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

                                 claqgb(N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in CGET04.

                           if ( EQUIL ) {
                              ROLDO = RCONDO
                              ROLDI = RCONDI
                           }

                           // Compute the 1-norm and infinity-norm of A.

                           ANORMO = CLANGB( '1', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )                            ANORMI = CLANGB( 'I', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )

                           // Factor the matrix A.

                           cgbtrf(N, N, KL, KU, AFB, LDAFB, IWORK, INFO );

                           // Form the inverse of A.

                           claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), WORK, LDB );
                           SRNAMT = 'CGBTRS'
                           cgbtrs('No transpose', N, KL, KU, N, AFB, LDAFB, IWORK, WORK, LDB, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = CLANGE( '1', N, N, WORK, LDB, RWORK )
                           if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDO = ONE
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM
                           }

                           // Compute the infinity-norm condition number
                           // of A.

                           AINVNM = CLANGE( 'I', N, N, WORK, LDB, RWORK )
                           if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDI = ONE
                           } else {
                              RCONDI = ( ONE / ANORMI ) / AINVNM
                           }
                        }

                        for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 90

                           // Do for each value of TRANS.

                           TRANS = TRANSS( ITRAN )
                           if ( ITRAN.EQ.1 ) {
                              RCONDC = RCONDO
                           } else {
                              RCONDC = RCONDI
                           }

                           // Restore the matrix A.

                           clacpy('Full', KL+KU+1, N, ASAV, LDA, A, LDA );

                           // Form an exact solution and set the right hand
                           // side.

                           SRNAMT = 'CLARHS'
                           clarhs(PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO );
                           XTYPE = 'C'
                           clacpy('Full', N, NRHS, B, LDB, BSAV, LDB );

                           if ( NOFACT .AND. ITRAN.EQ.1 ) {

                              // --- Test CGBSV  ---

                              // Compute the LU factorization of the matrix
                              // and solve the system.

                              clacpy('Full', KL+KU+1, N, A, LDA, AFB( KL+1 ), LDAFB )                               CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDB );

                              SRNAMT = 'CGBSV '
                              cgbsv(N, KL, KU, NRHS, AFB, LDAFB, IWORK, X, LDB, INFO );

                              // Check error code from CGBSV .

                              if (INFO.NE.IZERO) CALL ALAERH( PATH, 'CGBSV ', INFO, IZERO, ' ', N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                              // Reconstruct matrix from factors and
                              // compute residual.

                              cgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              NT = 1
                              if ( IZERO.EQ.0 ) {

                                 // Compute residual of the computed
                                 // solution.

                                 clacpy('Full', N, NRHS, B, LDB, WORK, LDB )                                  CALL CGBT02( 'No transpose', N, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

                                 // Check solution from generated exact
                                 // solution.

                                 cget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );
                                 NT = 3
                              }

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (K = 1; K <= NT; K++) { // 50
                                 if ( RESULT( K ).GE.THRESH ) {
                                    if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9997 )'CGBSV ', N, KL, KU, IMAT, K, RESULT( K );
                                    NFAIL = NFAIL + 1
                                 }
                              } // 50
                              NRUN = NRUN + NT
                           }

                           // --- Test CGBSVX ---

                           if (.NOT.PREFAC) CALL CLASET( 'Full', 2*KL+KU+1, N, CMPLX( ZERO ), CMPLX( ZERO ), AFB, LDAFB );
                           claset('Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDB );
                           if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                              // Equilibrate the matrix if FACT = 'F' and
                              // EQUED = 'R', 'C', or 'B'.

                              claqgb(N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                           }

                           // Solve the system and compute the condition
                           // number and error bounds using CGBSVX.

                           SRNAMT = 'CGBSVX'
                           cgbsvx(FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( LDB+1 ), B, LDB, X, LDB, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                           // Check the error code from CGBSVX.

                           if (INFO.NE.IZERO) CALL ALAERH( PATH, 'CGBSVX', INFO, IZERO, FACT // TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                           // Compare RWORK(2*NRHS+1) from CGBSVX with the
                           // computed reciprocal pivot growth RPVGRW

                           if ( INFO.NE.0 ) {
                              ANRMPV = ZERO
                              for (J = 1; J <= INFO; J++) { // 70
                                 DO 60 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )                                     ANRMPV = MAX( ANRMPV, ABS( A( I+( J-1 )*LDA ) ) )
                                 } // 60
                              } // 70
                              RPVGRW = CLANTB( 'M', 'U', 'N', INFO, MIN( INFO-1, KL+KU ), AFB( MAX( 1, KL+KU+2-INFO ) ), LDAFB, RDUM )
                              if ( RPVGRW.EQ.ZERO ) {
                                 RPVGRW = ONE
                              } else {
                                 RPVGRW = ANRMPV / RPVGRW
                              }
                           } else {
                              RPVGRW = CLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, RDUM )
                              if ( RPVGRW.EQ.ZERO ) {
                                 RPVGRW = ONE
                              } else {
                                 RPVGRW = CLANGB( 'M', N, KL, KU, A, LDA, RDUM ) / RPVGRW
                              }
                           }
                           RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) ) / MAX( RWORK( 2*NRHS+1 ), RPVGRW ) / SLAMCH( 'E' )

                           if ( .NOT.PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              cgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) );
                              K1 = 1
                           } else {
                              K1 = 2
                           }

                           if ( INFO.EQ.0 ) {
                              TRFCON = .FALSE.

                              // Compute residual of the computed solution.

                              clacpy('Full', N, NRHS, BSAV, LDB, WORK, LDB )                               CALL CGBT02( TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                              // Check solution from generated exact
                              // solution.

                              IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                                  CALL CGET04( N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) )
                              } else {
                                 if ( ITRAN.EQ.1 ) {
                                    ROLDC = ROLDO
                                 } else {
                                    ROLDC = ROLDI
                                 }
                                 cget04(N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) );
                              }

                              // Check the error bounds from iterative
                              // refinement.

                              cgbt05(TRANS, N, KL, KU, NRHS, ASAV, LDA, BSAV, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                           } else {
                              TRFCON = .TRUE.
                           }

                           // Compare RCOND from CGBSVX with the computed
                           // value in RCONDC.

                           RESULT( 6 ) = SGET06( RCOND, RCONDC )

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( .NOT.TRFCON ) {
                              for (K = K1; K <= NTESTS; K++) { // 80
                                 if ( RESULT( K ).GE.THRESH ) {
                                    if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                                    if ( PREFAC ) {
                                       WRITE( NOUT, FMT = 9995 ) 'CGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K )
                                    } else {
                                       WRITE( NOUT, FMT = 9996 ) 'CGBSVX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K )
                                    }
                                    NFAIL = NFAIL + 1
                                 }
                              } // 80
                              NRUN = NRUN + 7 - K1
                           } else {
                              IF( RESULT( 1 ).GE.THRESH .AND. .NOT. PREFAC ) THEN                                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                              if ( RESULT( 6 ).GE.THRESH ) {
                                 if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                              if ( RESULT( 7 ).GE.THRESH ) {
                                 if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                           }

                     // --- Test CGBSVXX ---

                     // Restore the matrices A and B.

                      // write(*,*) 'begin cgbsvxx testing'
                      clacpy('Full', KL+KU+1, N, ASAV, LDA, A, LDA );
                     clacpy('Full', N, NRHS, BSAV, LDB, B, LDB );
                      if (.NOT.PREFAC) CALL CLASET( 'Full', 2*KL+KU+1, N, CMPLX( ZERO ), CMPLX( ZERO ), AFB, LDAFB );
                     claset('Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDB );
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        claqgb(N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED );
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using CGBSVXX.

                     SRNAMT = 'CGBSVXX'
                     N_ERR_BNDS = 3
                     cgbsvxx(FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB, X, LDB, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, RWORK, INFO );

                     // Check the error code from CGBSVXX.

                     if (INFO.EQ.N+1) GOTO 90;
                     if ( INFO.NE.IZERO ) {
                        alaerh(PATH, 'CGBSVXX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GOTO 90
                     }

                     // Compare rpvgrw_svxx from CGESVXX with the computed
                     // reciprocal pivot growth factor RPVGRW


                     if ( INFO .GT. 0 .AND. INFO .LT. N+1 ) {
                        RPVGRW = CLA_GBRPVGRW(N, KL, KU, INFO, A, LDA, AFB, LDAFB)
                     } else {
                        RPVGRW = CLA_GBRPVGRW(N, KL, KU, N, A, LDA, AFB, LDAFB)
                     }
                      RESULT( 7 ) = ABS( RPVGRW-rpvgrw_svxx ) / MAX( rpvgrw_svxx, RPVGRW ) / SLAMCH( 'E' )

                     if ( .NOT.PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        cgbt01(N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     if ( INFO.EQ.0 ) {
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        clacpy('Full', N, NRHS, BSAV, LDB, WORK, LDB )                         CALL CGBT02( TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL CGET04( N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) )
                        } else {
                           if ( ITRAN.EQ.1 ) {
                              ROLDC = ROLDO
                           } else {
                              ROLDC = ROLDI
                           }
                           cget04(N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) );
                        }
                     } else {
                        TRFCON = .TRUE.
                     }

                     // Compare RCOND from CGBSVXX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = SGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( .NOT.TRFCON ) {
                        for (K = K1; K <= NTESTS; K++) { // 45
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
                        } // 45
                        NRUN = NRUN + 7 - K1
                     } else {
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 6 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 7 ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH );
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
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


      // Test Error Bounds from CGBSVXX

      cebchvxx(THRESH, PATH);

 9999 FORMAT( ' *** In CDRVGB, LA=', I5, ' is too small for N=', I5, ', KU=', I5, ', KL=', I5, / ' ==> Increase LA to at least ', I5 )
 9998 FORMAT( ' *** In CDRVGB, LAFB=', I5, ' is too small for N=', I5, ', KU=', I5, ', KL=', I5, / ' ==> Increase LAFB to at least ', I5 )
 9997 FORMAT( 1X, A, ', N=', I5, ', KL=', I5, ', KU=', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', I5, ',...), type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',', I5, ',...), EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=', G12.5 )

      RETURN

      // End of CDRVGBX

      }
