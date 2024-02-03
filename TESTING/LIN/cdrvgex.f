      SUBROUTINE CDRVGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      int                NTYPES;
      PARAMETER          ( NTYPES = 11 )
      int                NTESTS;
      PARAMETER          ( NTESTS = 7 )
      int                NTRAN;
      PARAMETER          ( NTRAN = 3 )
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, TRFCON, ZEROT;
      String             DIST, EQUED, FACT, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, ITRAN, IZERO, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT, N_ERR_BNDS;
      REAL               AINVNM, AMAX, ANORM, ANORMI, ANORMO, CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO, ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW, RPVGRW_SVXX
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANTR, SGET06, SLAMCH, CLA_GERPVGRW
      // EXTERNAL LSAME, CLANGE, CLANTR, SGET06, SLAMCH, CLA_GERPVGRW
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CERRVX, CGEEQU, CGESV, CGESVX, CGET01, CGET02, CGET04, CGET07, CGETRF, CGETRI, CLACPY, CLAQGE, CLARHS, CLASET, CLATB4, CLATMS, XLAENV, CGESVXX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX
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

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL CERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 90 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 80 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
            IF( ZEROT .AND. N.LT.IMAT-4 ) GO TO 80

            // Set up parameters with CLATB4 and generate a test matrix
            // with CLATMS.

            CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
            RCONDC = ONE / CNDNUM

            SRNAMT = 'CLATMS'
            CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO )

            // Check error code from CLATMS.

            IF( INFO.NE.0 ) THEN
               CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 80
            END IF

            // For types 5-7, zero one or more columns of the matrix to
           t // est that INFO is returned correctly.

            IF( ZEROT ) THEN
               IF( IMAT.EQ.5 ) THEN
                  IZERO = 1
               ELSE IF( IMAT.EQ.6 ) THEN
                  IZERO = N
               ELSE
                  IZERO = N / 2 + 1
               END IF
               IOFF = ( IZERO-1 )*LDA
               IF( IMAT.LT.7 ) THEN
                  DO 20 I = 1, N
                     A( IOFF+I ) = ZERO
   20             CONTINUE
               ELSE
                  CALL CLASET( 'Full', N, N-IZERO+1, CMPLX( ZERO ), CMPLX( ZERO ), A( IOFF+1 ), LDA )
               END IF
            ELSE
               IZERO = 0
            END IF

            // Save a copy of the matrix A in ASAV.

            CALL CLACPY( 'Full', N, N, A, LDA, ASAV, LDA )

            DO 70 IEQUED = 1, 4
               EQUED = EQUEDS( IEQUED )
               IF( IEQUED.EQ.1 ) THEN
                  NFACT = 3
               ELSE
                  NFACT = 1
               END IF

               DO 60 IFACT = 1, NFACT
                  FACT = FACTS( IFACT )
                  PREFAC = LSAME( FACT, 'F' )
                  NOFACT = LSAME( FACT, 'N' )
                  EQUIL = LSAME( FACT, 'E' )

                  IF( ZEROT ) THEN
                     IF( PREFAC ) GO TO 60
                     RCONDO = ZERO
                     RCONDI = ZERO

                  ELSE IF( .NOT.NOFACT ) THEN

                     // Compute the condition number for comparison with
                    t // he value returned by CGESVX (FACT = 'N' reuses
                    t // he condition number from the previous iteration
                     // with FACT = 'F').

                     CALL CLACPY( 'Full', N, N, ASAV, LDA, AFAC, LDA )
                     IF( EQUIL .OR. IEQUED.GT.1 ) THEN

                        // Compute row and column scale factors to
                        // equilibrate the matrix A.

                        CALL CGEEQU( N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, INFO )
                        IF( INFO.EQ.0 .AND. N.GT.0 ) THEN
                           IF( LSAME( EQUED, 'R' ) ) THEN
                              ROWCND = ZERO
                              COLCND = ONE
                           ELSE IF( LSAME( EQUED, 'C' ) ) THEN
                              ROWCND = ONE
                              COLCND = ZERO
                           ELSE IF( LSAME( EQUED, 'B' ) ) THEN
                              ROWCND = ZERO
                              COLCND = ZERO
                           END IF

                           // Equilibrate the matrix.

                           CALL CLAQGE( N, N, AFAC, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                        END IF
                     END IF

                     // Save the condition number of the non-equilibrated
                     // system for use in CGET04.

                     IF( EQUIL ) THEN
                        ROLDO = RCONDO
                        ROLDI = RCONDI
                     END IF

                     // Compute the 1-norm and infinity-norm of A.

                     ANORMO = CLANGE( '1', N, N, AFAC, LDA, RWORK )
                     ANORMI = CLANGE( 'I', N, N, AFAC, LDA, RWORK )

                     // Factor the matrix A.

                     CALL CGETRF( N, N, AFAC, LDA, IWORK, INFO )

                     // Form the inverse of A.

                     CALL CLACPY( 'Full', N, N, AFAC, LDA, A, LDA )
                     LWORK = NMAX*MAX( 3, NRHS )
                     CALL CGETRI( N, A, LDA, IWORK, WORK, LWORK, INFO )

                     // Compute the 1-norm condition number of A.

                     AINVNM = CLANGE( '1', N, N, A, LDA, RWORK )
                     IF( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDO = ONE
                     ELSE
                        RCONDO = ( ONE / ANORMO ) / AINVNM
                     END IF

                     // Compute the infinity-norm condition number of A.

                     AINVNM = CLANGE( 'I', N, N, A, LDA, RWORK )
                     IF( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDI = ONE
                     ELSE
                        RCONDI = ( ONE / ANORMI ) / AINVNM
                     END IF
                  END IF

                  DO 50 ITRAN = 1, NTRAN

                     // Do for each value of TRANS.

                     TRANS = TRANSS( ITRAN )
                     IF( ITRAN.EQ.1 ) THEN
                        RCONDC = RCONDO
                     ELSE
                        RCONDC = RCONDI
                     END IF

                     // Restore the matrix A.

                     CALL CLACPY( 'Full', N, N, ASAV, LDA, A, LDA )

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'CLARHS'
                     CALL CLARHS( PATH, XTYPE, 'Full', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     XTYPE = 'C'
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )

                     IF( NOFACT .AND. ITRAN.EQ.1 ) THEN

                        // --- Test CGESV  ---

                        // Compute the LU factorization of the matrix and
                        // solve the system.

                        CALL CLACPY( 'Full', N, N, A, LDA, AFAC, LDA )
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                        SRNAMT = 'CGESV '
                        CALL CGESV( N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO )

                        // Check error code from CGESV .

                        IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'CGESV ', INFO, IZERO, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL CGET01( N, N, A, LDA, AFAC, LDA, IWORK, RWORK, RESULT( 1 ) )
                        NT = 1
                        IF( IZERO.EQ.0 ) THEN

                           // Compute residual of the computed solution.

                           CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                            CALL CGET02( 'No transpose', N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

                           // Check solution from generated exact solution.

                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                           NT = 3
                        END IF

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 30 K = 1, NT
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'CGESV ', N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   30                   CONTINUE
                        NRUN = NRUN + NT
                     END IF

                     // --- Test CGESVX ---

                     IF( .NOT.PREFAC ) CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ZERO ), AFAC, LDA )
                     CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDA )
                     IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        CALL CLAQGE( N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                     END IF

                     // Solve the system and compute the condition number
                     // and error bounds using CGESVX.

                     SRNAMT = 'CGESVX'
                     CALL CGESVX( FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                     // Check the error code from CGESVX.

                     IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'CGESVX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                     // Compare RWORK(2*NRHS+1) from CGESVX with the
                     // computed reciprocal pivot growth factor RPVGRW

                     IF( INFO.NE.0 ) THEN
                        RPVGRW = CLANTR( 'M', 'U', 'N', INFO, INFO, AFAC, LDA, RDUM )
                        IF( RPVGRW.EQ.ZERO ) THEN
                           RPVGRW = ONE
                        ELSE
                           RPVGRW = CLANGE( 'M', N, INFO, A, LDA, RDUM ) / RPVGRW
                        END IF
                     ELSE
                        RPVGRW = CLANTR( 'M', 'U', 'N', N, N, AFAC, LDA, RDUM )
                        IF( RPVGRW.EQ.ZERO ) THEN
                           RPVGRW = ONE
                        ELSE
                           RPVGRW = CLANGE( 'M', N, N, A, LDA, RDUM ) / RPVGRW
                        END IF
                     END IF
                     RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) ) / MAX( RWORK( 2*NRHS+1 ), RPVGRW ) / SLAMCH( 'E' )

                     IF( .NOT.PREFAC ) THEN

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL CGET01( N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF

                     IF( INFO.EQ.0 ) THEN
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL CGET02( TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        ELSE
                           IF( ITRAN.EQ.1 ) THEN
                              ROLDC = ROLDO
                           ELSE
                              ROLDC = ROLDI
                           END IF
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) )
                        END IF

                        // Check the error bounds from iterative
                        // refinement.

                        CALL CGET07( TRANS, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, .TRUE., RWORK( NRHS+1 ), RESULT( 4 ) )
                     ELSE
                        TRFCON = .TRUE.
                     END IF

                     // Compare RCOND from CGESVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = SGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     IF( .NOT.TRFCON ) THEN
                        DO 40 K = K1, NTESTS
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9997 )'CGESVX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K )
                              ELSE
                                 WRITE( NOUT, FMT = 9998 )'CGESVX', FACT, TRANS, N, IMAT, K, RESULT( K )
                              END IF
                              NFAIL = NFAIL + 1
                           END IF
   40                   CONTINUE
                        NRUN = NRUN + 7 - K1
                     ELSE
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVX', FACT, TRANS, N, IMAT, 1, RESULT( 1 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 6 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVX', FACT, TRANS, N, IMAT, 6, RESULT( 6 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 7 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVX', FACT, TRANS, N, IMAT, 7, RESULT( 7 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF

                     END IF

                     // --- Test CGESVXX ---

                     // Restore the matrices A and B.


                     CALL CLACPY( 'Full', N, N, ASAV, LDA, A, LDA )
                     CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, B, LDA )
                      IF( .NOT.PREFAC ) CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ZERO ), AFAC, LDA )
                     CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDA )
                     IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN

                        // Equilibrate the matrix if FACT = 'F' and
                        // EQUED = 'R', 'C', or 'B'.

                        CALL CLAQGE( N, N, A, LDA, S, S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                     END IF

                     // Solve the system and compute the condition number
                     // and error bounds using CGESVXX.

                     SRNAMT = 'CGESVXX'
                     N_ERR_BNDS = 3
                     CALL CGESVXX( FACT, TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, S, S( N+1 ), B, LDA, X, LDA, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, RWORK, INFO )

                     // Check the error code from CGESVXX.

                     IF( INFO.EQ.N+1 ) GOTO 50
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'CGESVXX', INFO, IZERO, FACT // TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GOTO 50
                     END IF

                     // Compare rpvgrw_svxx from CGESVXX with the computed
                     // reciprocal pivot growth factor RPVGRW


                     IF ( INFO .GT. 0 .AND. INFO .LT. N+1 ) THEN
                        RPVGRW = CLA_GERPVGRW (N, INFO, A, LDA, AFAC, LDA)
                     ELSE
                        RPVGRW = CLA_GERPVGRW (N, N, A, LDA, AFAC, LDA)
                     ENDIF
                      RESULT( 7 ) = ABS( RPVGRW-rpvgrw_svxx ) / MAX( rpvgrw_svxx, RPVGRW ) / SLAMCH( 'E' )

                     IF( .NOT.PREFAC ) THEN

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL CGET01( N, N, A, LDA, AFAC, LDA, IWORK, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF

                     IF( INFO.EQ.0 ) THEN
                        TRFCON = .FALSE.

                        // Compute residual of the computed solution.

                        CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL CGET02( TRANS, N, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        ELSE
                           IF( ITRAN.EQ.1 ) THEN
                              ROLDC = ROLDO
                           ELSE
                              ROLDC = ROLDI
                           END IF
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) )
                        END IF
                     ELSE
                        TRFCON = .TRUE.
                     END IF

                     // Compare RCOND from CGESVXX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = SGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     IF( .NOT.TRFCON ) THEN
                        DO 45 K = K1, NTESTS
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9997 )'CGESVXX', FACT, TRANS, N, EQUED, IMAT, K, RESULT( K )
                              ELSE
                                 WRITE( NOUT, FMT = 9998 )'CGESVXX', FACT, TRANS, N, IMAT, K, RESULT( K )
                              END IF
                              NFAIL = NFAIL + 1
                           END IF
 45                     CONTINUE
                        NRUN = NRUN + 7 - K1
                     ELSE
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC ) THEN                            IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVXX', FACT, TRANS, N, EQUED, IMAT, 1, RESULT( 1 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVXX', FACT, TRANS, N, IMAT, 1, RESULT( 1 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 6 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVXX', FACT, TRANS, N, EQUED, IMAT, 6, RESULT( 6 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVXX', FACT, TRANS, N, IMAT, 6, RESULT( 6 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 7 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CGESVXX', FACT, TRANS, N, EQUED, IMAT, 7, RESULT( 7 )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CGESVXX', FACT, TRANS, N, IMAT, 7, RESULT( 7 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF

                     END IF

   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )


      // Test Error Bounds for CGESVXX

      CALL CEBCHVXX(THRESH, PATH)

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test(', I2, ') =',
     $      G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5,
     $      ', type ', I2, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5,
     $      ', EQUED=''', A1, ''', type ', I2, ', test(', I1, ')=',
     $      G12.5 )
      RETURN

      // End of CDRVGEX

      END
