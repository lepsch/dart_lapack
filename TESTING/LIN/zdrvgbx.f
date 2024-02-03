      SUBROUTINE ZDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA, AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

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
      double             RWORK( * ), S( * );
      COMPLEX*16         A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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
      double             AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, RPVGRW_SVXX;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLAMCH, ZLANGB, ZLANGE, ZLANTB, ZLA_GBRPVGRW       EXTERNAL           LSAME, DGET06, DLAMCH, ZLANGB, ZLANGE, ZLANTB, ZLA_GBRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGBEQU, ZGBSV, ZGBSVX, ZGBT01, ZGBT02, ZGBT05, ZGBTRF, ZGBTRS, ZGET04, ZLACPY, ZLAQGB, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZGBSVXX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, MAX, MIN
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
      PATH( 2: 3 ) = 'GB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDB = MAX( N, 1 )
         XTYPE = 'N'

         // Set limits on the number of loop iterations.

         NKL = MAX( 1, MIN( N, 4 ) )
         IF( N.EQ.0 ) NKL = 1
         NKU = NKL
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 140 IKL = 1, NKL

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
            DO 130 IKU = 1, NKU

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
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
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

               DO 120 IMAT = 1, NIMAT

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

                  // Skip types 2, 3, or 4 if the matrix is too small.

                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 ) GO TO 120

                  // Set up parameters with ZLATB4 and generate a
                 t // est matrix with ZLATMS.

                  CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
                  RCONDC = ONE / CNDNUM

                  SRNAMT = 'ZLATMS'
                  CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK, INFO )

                  // Check the error code from ZLATMS.

                  if ( INFO.NE.0 ) {
                     CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 120
                  }

                  // For types 2, 3, and 4, zero one or more columns of
                 t // he matrix to test that INFO is returned correctly.

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
                        DO 20 I = I1, I2
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                     } else {
                        DO 40 J = IZERO, N
                           DO 30 I = MAX( 1, KU+2-J ), MIN( KL+KU+1, KU+1+( N-J ) )
                              A( IOFF+I ) = ZERO
   30                      CONTINUE
                           IOFF = IOFF + LDA
   40                   CONTINUE
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  CALL ZLACPY( 'Full', KL+KU+1, N, A, LDA, ASAV, LDA )

                  DO 110 IEQUED = 1, 4
                     EQUED = EQUEDS( IEQUED )
                     if ( IEQUED.EQ.1 ) {
                        NFACT = 3
                     } else {
                        NFACT = 1
                     }

                     DO 100 IFACT = 1, NFACT
                        FACT = FACTS( IFACT )
                        PREFAC = LSAME( FACT, 'F' )
                        NOFACT = LSAME( FACT, 'N' )
                        EQUIL = LSAME( FACT, 'E' )

                        if ( ZEROT ) {
                           IF( PREFAC ) GO TO 100
                           RCONDO = ZERO
                           RCONDI = ZERO

                        } else if ( .NOT.NOFACT ) {

                           // Compute the condition number for comparison
                           // with the value returned by DGESVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           CALL ZLACPY( 'Full', KL+KU+1, N, ASAV, LDA, AFB( KL+1 ), LDAFB )
                           if ( EQUIL .OR. IEQUED.GT.1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              CALL ZGBEQU( N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO )
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

                                 CALL ZLAQGB( N, N, KL, KU, AFB( KL+1 ), LDAFB, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in ZGET04.

                           if ( EQUIL ) {
                              ROLDO = RCONDO
                              ROLDI = RCONDI
                           }

                           // Compute the 1-norm and infinity-norm of A.

                           ANORMO = ZLANGB( '1', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )                            ANORMI = ZLANGB( 'I', N, KL, KU, AFB( KL+1 ), LDAFB, RWORK )

                           // Factor the matrix A.

                           CALL ZGBTRF( N, N, KL, KU, AFB, LDAFB, IWORK, INFO )

                           // Form the inverse of A.

                           CALL ZLASET( 'Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, LDB )
                           SRNAMT = 'ZGBTRS'
                           CALL ZGBTRS( 'No transpose', N, KL, KU, N, AFB, LDAFB, IWORK, WORK, LDB, INFO )

                           // Compute the 1-norm condition number of A.

                           AINVNM = ZLANGE( '1', N, N, WORK, LDB, RWORK )
                           if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDO = ONE
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM
                           }

                           // Compute the infinity-norm condition number
                           // of A.

                           AINVNM = ZLANGE( 'I', N, N, WORK, LDB, RWORK )
                           if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDI = ONE
                           } else {
                              RCONDI = ( ONE / ANORMI ) / AINVNM
                           }
                        }

                        DO 90 ITRAN = 1, NTRAN

                           // Do for each value of TRANS.

                           TRANS = TRANSS( ITRAN )
                           if ( ITRAN.EQ.1 ) {
                              RCONDC = RCONDO
                           } else {
                              RCONDC = RCONDI
                           }

                           // Restore the matrix A.

                           CALL ZLACPY( 'Full', KL+KU+1, N, ASAV, LDA, A, LDA )

                           // Form an exact solution and set the right hand
                           // side.

                           SRNAMT = 'ZLARHS'
                           CALL ZLARHS( PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO )
                           XTYPE = 'C'
                           CALL ZLACPY( 'Full', N, NRHS, B, LDB, BSAV, LDB )

                           if ( NOFACT .AND. ITRAN.EQ.1 ) {

                              // --- Test ZGBSV  ---

                              // Compute the LU factorization of the matrix
                              // and solve the system.

                              CALL ZLACPY( 'Full', KL+KU+1, N, A, LDA, AFB( KL+1 ), LDAFB )                               CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDB )

                              SRNAMT = 'ZGBSV '
                              CALL ZGBSV( N, KL, KU, NRHS, AFB, LDAFB, IWORK, X, LDB, INFO )

                              // Check error code from ZGBSV .

                              IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'ZGBSV ', INFO, IZERO, ' ', N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT )

                              // Reconstruct matrix from factors and
                              // compute residual.

                              CALL ZGBT01( N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) )
                              NT = 1
                              if ( IZERO.EQ.0 ) {

                                 // Compute residual of the computed
                                 // solution.

                                 CALL ZLACPY( 'Full', N, NRHS, B, LDB, WORK, LDB )                                  CALL ZGBT02( 'No transpose', N, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) )

                                 // Check solution from generated exact
                                 // solution.

                                 CALL ZGET04( N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) )
                                 NT = 3
                              }

                              // Print information about the tests that did
                              // not pass the threshold.

                              DO 50 K = 1, NT
                                 if ( RESULT( K ).GE.THRESH ) {
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9997 )'ZGBSV ', N, KL, KU, IMAT, K, RESULT( K )
                                    NFAIL = NFAIL + 1
                                 }
   50                         CONTINUE
                              NRUN = NRUN + NT
                           }

                           // --- Test ZGBSVX ---

                           IF( .NOT.PREFAC ) CALL ZLASET( 'Full', 2*KL+KU+1, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFB, LDAFB )
                           CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDB )
                           if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                              // Equilibrate the matrix if FACT = 'F' and
                              // EQUED = 'R', 'C', or 'B'.

                              CALL ZLAQGB( N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                           }

                           // Solve the system and compute the condition
                           // number and error bounds using ZGBSVX.

                           SRNAMT = 'ZGBSVX'
                           CALL ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( LDB+1 ), B, LDB, X, LDB, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                           // Check the error code from ZGBSVX.

                           IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'ZGBSVX', INFO, IZERO, FACT // TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT )

                           // Compare RWORK(2*NRHS+1) from ZGBSVX with the
                           // computed reciprocal pivot growth RPVGRW

                           if ( INFO.NE.0 ) {
                              ANRMPV = ZERO
                              DO 70 J = 1, INFO
                                 DO 60 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )                                     ANRMPV = MAX( ANRMPV, ABS( A( I+( J-1 )*LDA ) ) )
   60                            CONTINUE
   70                         CONTINUE
                              RPVGRW = ZLANTB( 'M', 'U', 'N', INFO, MIN( INFO-1, KL+KU ), AFB( MAX( 1, KL+KU+2-INFO ) ), LDAFB, RDUM )
                              if ( RPVGRW.EQ.ZERO ) {
                                 RPVGRW = ONE
                              } else {
                                 RPVGRW = ANRMPV / RPVGRW
                              }
                           } else {
                              RPVGRW = ZLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, RDUM )
                              if ( RPVGRW.EQ.ZERO ) {
                                 RPVGRW = ONE
                              } else {
                                 RPVGRW = ZLANGB( 'M', N, KL, KU, A, LDA, RDUM ) / RPVGRW
                              }
                           }
                           RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) ) / MAX( RWORK( 2*NRHS+1 ), RPVGRW ) / DLAMCH( 'E' )

                           if ( .NOT.PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              CALL ZGBT01( N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK, RESULT( 1 ) )
                              K1 = 1
                           } else {
                              K1 = 2
                           }

                           if ( INFO.EQ.0 ) {
                              TRFCON = .FALSE.

                              // Compute residual of the computed solution.

                              CALL ZLACPY( 'Full', N, NRHS, BSAV, LDB, WORK, LDB )                               CALL ZGBT02( TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                              // Check solution from generated exact
                              // solution.

                              IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                                  CALL ZGET04( N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) )
                              } else {
                                 if ( ITRAN.EQ.1 ) {
                                    ROLDC = ROLDO
                                 } else {
                                    ROLDC = ROLDI
                                 }
                                 CALL ZGET04( N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) )
                              }

                              // Check the error bounds from iterative
                              // refinement.

                              CALL ZGBT05( TRANS, N, KL, KU, NRHS, ASAV, LDA, BSAV, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                           } else {
                              TRFCON = .TRUE.
                           }

                           // Compare RCOND from ZGBSVX with the computed
                           // value in RCONDC.

                           RESULT( 6 ) = DGET06( RCOND, RCONDC )

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( .NOT.TRFCON ) {
                              DO 80 K = K1, NTESTS
                                 if ( RESULT( K ).GE.THRESH ) {
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                                    if ( PREFAC ) {
                                       WRITE( NOUT, FMT = 9995 ) 'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K )
                                    } else {
                                       WRITE( NOUT, FMT = 9996 ) 'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K )
                                    }
                                    NFAIL = NFAIL + 1
                                 }
   80                         CONTINUE
                              NRUN = NRUN + 7 - K1
                           } else {
                              IF( RESULT( 1 ).GE.THRESH .AND. .NOT. PREFAC ) THEN                                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                              if ( RESULT( 6 ).GE.THRESH ) {
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                              if ( RESULT( 7 ).GE.THRESH ) {
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                                 if ( PREFAC ) {
                                    WRITE( NOUT, FMT = 9995 )'ZGBSVX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 )
                                 } else {
                                    WRITE( NOUT, FMT = 9996 )'ZGBSVX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 )
                                 }
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              }
                           }

                     // --- Test ZGBSVXX ---

                     // Restore the matrices A and B.

                      // write(*,*) 'begin zgbsvxx testing'
                      CALL ZLACPY( 'Full', KL+KU+1, N, ASAV, LDA, A, LDA )
                     CALL ZLACPY( 'Full', N, NRHS, BSAV, LDB, B, LDB )
                      IF( .NOT.PREFAC ) CALL ZLASET( 'Full', 2*KL+KU+1, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFB, LDAFB )
                     CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDB )
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        CALL ZLAQGB( N, N, KL, KU, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZGBSVXX.

                     SRNAMT = 'ZGBSVXX'
                     N_ERR_BNDS = 3
                     CALL ZGBSVXX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB, X, LDB, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, RWORK, INFO )

                     // Check the error code from ZGBSVXX.

                     IF( INFO.EQ.N+1 ) GOTO 90
                     if ( INFO.NE.IZERO ) {
                        CALL ALAERH( PATH, 'ZGBSVXX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GOTO 90
                     }

                     // Compare rpvgrw_svxx from ZGESVXX with the computed
                     // reciprocal pivot growth factor RPVGRW


                     if ( INFO .GT. 0 .AND. INFO .LT. N+1 ) {
                        RPVGRW = ZLA_GBRPVGRW(N, KL, KU, INFO, A, LDA, AFB, LDAFB)
                     } else {
                        RPVGRW = ZLA_GBRPVGRW(N, KL, KU, N, A, LDA, AFB, LDAFB)
                     ENDIF
                      RESULT( 7 ) = ABS( RPVGRW-rpvgrw_svxx ) / MAX( rpvgrw_svxx, RPVGRW ) / DLAMCH( 'E' )

                     if ( .NOT.PREFAC ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL ZGBT01( N, N, KL, KU, A, LDA, AFB, LDAFB, IWORK, WORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     if ( INFO.EQ.0 ) {
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        CALL ZLACPY( 'Full', N, NRHS, BSAV, LDB, WORK, LDB )                         CALL ZGBT02( TRANS, N, N, KL, KU, NRHS, ASAV, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL ZGET04( N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) )
                        } else {
                           if ( ITRAN.EQ.1 ) {
                              ROLDC = ROLDO
                           } else {
                              ROLDC = ROLDI
                           }
                           CALL ZGET04( N, NRHS, X, LDB, XACT, LDB, ROLDC, RESULT( 3 ) )
                        }
                     } else {
                        TRFCON = .TRUE.
                     }

                     // Compare RCOND from ZGBSVXX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     if ( .NOT.TRFCON ) {
                        DO 45 K = K1, NTESTS
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
 45                     CONTINUE
                        NRUN = NRUN + 7 - K1
                     } else {
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 1, RESULT( 1 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 1, RESULT( 1 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 6 ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 6, RESULT( 6 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 6, RESULT( 6 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }
                        if ( RESULT( 7 ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9995 )'ZGBSVXX', FACT, TRANS, N, KL, KU, EQUED, IMAT, 7, RESULT( 7 )
                           } else {
                              WRITE( NOUT, FMT = 9996 )'ZGBSVXX', FACT, TRANS, N, KL, KU, IMAT, 7, RESULT( 7 )
                           }
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        }

                     }

   90                   CONTINUE
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )


      // Test Error Bounds from ZGBSVXX

      CALL ZEBCHVXX(THRESH, PATH)

 9999 FORMAT( ' *** In ZDRVGB, LA=', I5, ' is too small for N=', I5,
     $      ', KU=', I5, ', KL=', I5, / ' ==> Increase LA to at least ',
     $      I5 )
 9998 FORMAT( ' *** In ZDRVGB, LAFB=', I5, ' is too small for N=', I5,
     $      ', KU=', I5, ', KL=', I5, /
     $      ' ==> Increase LAFB to at least ', I5 )
 9997 FORMAT( 1X, A, ', N=', I5, ', KL=', I5, ', KU=', I5, ', type ',
     $      I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',',
     $      I5, ',...), type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',',
     $      I5, ',...), EQUED=''', A1, ''', type ', I1, ', test(', I1,
     $      ')=', G12.5 )

      RETURN

      // End of ZDRVGBX

      }
