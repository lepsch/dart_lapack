      SUBROUTINE ZDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      int                NTYPES, NTESTS;
      PARAMETER          ( NTYPES = 11, NTESTS = 6 )
      int                NFACT;
      PARAMETER          ( NFACT = 2 )
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, EQUED, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT, N_ERR_BNDS;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC, RPVGRW_SVXX;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANSY;
      // EXTERNAL DGET06, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGET04, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZLATSY, ZPOT05, ZSYSV, ZSYSVX, ZSYT01, ZSYT02, ZSYTRF, ZSYTRI2, ZSYSVXX
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
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..
*
      // Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'SY'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      LWORK = MAX( 2*NMAX, NMAX*NRHS )
*
      // Test the error exits
*
      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0
*
      // Set the block size and minimum block size for testing.
*
      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )
*
      // Do for each value of N in NVAL
*
      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1
*
         DO 170 IMAT = 1, NIMAT
*
            // Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170
*
            // Skip types 3, 4, 5, or 6 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 170
*
            // Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 160 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
               IF( IMAT.NE.NTYPES ) THEN
*
                  // Set up parameters with ZLATB4 and generate a test
                  // matrix with ZLATMS.
*
                  CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
                  SRNAMT = 'ZLATMS'
                  CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )
*
                  // Check error code from ZLATMS.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 160
                  END IF
*
                  // For types 3-6, zero one or more rows and columns of
                 t // he matrix to test that INFO is returned correctly.
*
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.3 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.4 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
                     IF( IMAT.LT.6 ) THEN
*
                        // Set row and column IZERO to zero.
*
                        IF( IUPLO.EQ.1 ) THEN
                           IOFF = ( IZERO-1 )*LDA
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        ELSE
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        END IF
                     ELSE
                        IF( IUPLO.EQ.1 ) THEN
*
                           // Set the first IZERO rows to zero.
*
                           IOFF = 0
                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = ZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                        ELSE
*
                           // Set the last IZERO rows to zero.
*
                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = ZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        END IF
                     END IF
                  ELSE
                     IZERO = 0
                  END IF
               ELSE
*
                  // IMAT = NTYPES:  Use a special block diagonal matrix to
                 t // est alternate code for the 2-by-2 blocks.
*
                  CALL ZLATSY( UPLO, N, A, LDA, ISEED )
               END IF
*
               DO 150 IFACT = 1, NFACT
*
                  // Do first for FACT = 'F', then for other values.
*
                  FACT = FACTS( IFACT )
*
                  // Compute the condition number for comparison with
                 t // he value returned by ZSYSVX.
*
                  IF( ZEROT ) THEN
                     IF( IFACT.EQ.1 ) GO TO 150
                     RCONDC = ZERO
*
                  ELSE IF( IFACT.EQ.1 ) THEN
*
                     // Compute the 1-norm of A.
*
                     ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK )
*
                     // Factor the matrix A.
*
                     CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL ZSYTRF( UPLO, N, AFAC, LDA, IWORK, WORK, LWORK, INFO )
*
                     // Compute inv(A) and take its norm.
*
                     CALL ZLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     LWORK = (N+NB+1)*(NB+3)
                     CALL ZSYTRI2( UPLO, N, AINV, LDA, IWORK, WORK, LWORK, INFO )
                     AINVNM = ZLANSY( '1', UPLO, N, AINV, LDA, RWORK )
*
                     // Compute the 1-norm condition number of A.
*
                     IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
                  END IF
*
                  // Form an exact solution and set the right hand side.
*
                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                  XTYPE = 'C'
*
                  // --- Test ZSYSV  ---
*
                  IF( IFACT.EQ.2 ) THEN
                     CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     // Factor the matrix and solve the system using ZSYSV.
*
                     SRNAMT = 'ZSYSV '
                     CALL ZSYSV( UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO )
*
                     // Adjust the expected value of INFO to account for
                     // pivoting.
*
                     K = IZERO
                     IF( K.GT.0 ) THEN
  100                   CONTINUE
                        IF( IWORK( K ).LT.0 ) THEN
                           IF( IWORK( K ).NE.-K ) THEN
                              K = -IWORK( K )
                              GO TO 100
                           END IF
                        ELSE IF( IWORK( K ).NE.K ) THEN
                           K = IWORK( K )
                           GO TO 100
                        END IF
                     END IF
*
                     // Check error code from ZSYSV .
*
                     IF( INFO.NE.K ) THEN
                        CALL ALAERH( PATH, 'ZSYSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 120
                     ELSE IF( INFO.NE.0 ) THEN
                        GO TO 120
                     END IF
*
                     // Reconstruct matrix from factors and compute
                     // residual.
*
                     CALL ZSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) )
*
                     // Compute residual of the computed solution.
*
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )
*
                     // Check solution from generated exact solution.
*
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                     NT = 3
*
                     // Print information about the tests that did not pass
                    t // he threshold.
*
                     DO 110 K = 1, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'ZSYSV ', UPLO, N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  END IF
*
                  // --- Test ZSYSVX ---
*
                  IF( IFACT.EQ.2 ) CALL ZLASET( UPLO, N, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, LDA )
                  CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA )
*
                  // Solve the system and compute the condition number and
                  // error bounds using ZSYSVX.
*
                  SRNAMT = 'ZSYSVX'
                  CALL ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, LWORK, RWORK( 2*NRHS+1 ), INFO )
*
                  // Adjust the expected value of INFO to account for
                  // pivoting.
*
                  K = IZERO
                  IF( K.GT.0 ) THEN
  130                CONTINUE
                     IF( IWORK( K ).LT.0 ) THEN
                        IF( IWORK( K ).NE.-K ) THEN
                           K = -IWORK( K )
                           GO TO 130
                        END IF
                     ELSE IF( IWORK( K ).NE.K ) THEN
                        K = IWORK( K )
                        GO TO 130
                     END IF
                  END IF
*
                  // Check the error code from ZSYSVX.
*
                  IF( INFO.NE.K ) THEN
                     CALL ALAERH( PATH, 'ZSYSVX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 150
                  END IF
*
                  IF( INFO.EQ.0 ) THEN
                     IF( IFACT.GE.2 ) THEN
*
                        // Reconstruct matrix from factors and compute
                        // residual.
*
                        CALL ZSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF
*
                     // Compute residual of the computed solution.
*
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )
*
                     // Check solution from generated exact solution.
*
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
*
                     // Check the error bounds from iterative refinement.
*
                     CALL ZPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
*
                  // Compare RCOND from ZSYSVX with the computed value
                  // in RCONDC.
*
                  RESULT( 6 ) = DGET06( RCOND, RCONDC )
*
                  // Print information about the tests that did not pass
                 t // he threshold.
*
                  DO 140 K = K1, 6
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'ZSYSVX', FACT, UPLO, N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  140             CONTINUE
                  NRUN = NRUN + 7 - K1
*
                  // --- Test ZSYSVXX ---
*
                  // Restore the matrices A and B.
*
                  IF( IFACT.EQ.2 ) CALL ZLASET( UPLO, N, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, LDA )
                  CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA )
*
                  // Solve the system and compute the condition number
                  // and error bounds using ZSYSVXX.
*
                  SRNAMT = 'ZSYSVXX'
                  N_ERR_BNDS = 3
                  EQUED = 'N'
                  CALL ZSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, WORK( N+1 ), B, LDA, X, LDA, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, RWORK, INFO )
*
                  // Adjust the expected value of INFO to account for
                  // pivoting.
*
                  K = IZERO
                  IF( K.GT.0 ) THEN
 135                 CONTINUE
                     IF( IWORK( K ).LT.0 ) THEN
                        IF( IWORK( K ).NE.-K ) THEN
                           K = -IWORK( K )
                           GO TO 135
                        END IF
                     ELSE IF( IWORK( K ).NE.K ) THEN
                        K = IWORK( K )
                        GO TO 135
                     END IF
                  END IF
*
                  // Check the error code from ZSYSVXX.
*
                  IF( INFO.NE.K .AND. INFO.LE.N ) THEN
                     CALL ALAERH( PATH, 'ZSYSVXX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 150
                  END IF
*
                  IF( INFO.EQ.0 ) THEN
                     IF( IFACT.GE.2 ) THEN
*
                  // Reconstruct matrix from factors and compute
                  // residual.
*
                        CALL ZSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK(2*NRHS+1), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF
*
                  // Compute residual of the computed solution.
*
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )
                     RESULT( 2 ) = 0.0
*
                  // Check solution from generated exact solution.
*
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
*
                  // Check the error bounds from iterative refinement.
*
                     CALL ZPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
*
                  // Compare RCOND from ZSYSVXX with the computed value
                  // in RCONDC.
*
                  RESULT( 6 ) = DGET06( RCOND, RCONDC )
*
                  // Print information about the tests that did not pass
                 t // he threshold.
*
                  DO 85 K = K1, 6
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'ZSYSVXX', FACT, UPLO, N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
 85               CONTINUE
                  NRUN = NRUN + 7 - K1
*
  150          CONTINUE
*
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
*
      // Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*

      // Test Error Bounds from ZSYSVXX

      CALL ZEBCHVXX(THRESH, PATH)

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5,
     $      ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN
*
      // End of ZDRVSYX
*
      END
