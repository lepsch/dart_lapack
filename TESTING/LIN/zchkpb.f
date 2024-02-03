      SUBROUTINE ZCHKPB( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NBVAL( * ), NSVAL( * ), NVAL( * );
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
      PARAMETER          ( NTYPES = 8, NTESTS = 7 )
      int                NBW;
      PARAMETER          ( NBW = 4 )
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IKD, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IW, IZERO, K, KD, KL, KOFF, KU, LDA, LDAB, MODE, N, NB, NERRS, NFAIL, NIMAT, NKD, NRHS, NRUN;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), KDVAL( NBW );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANGE, ZLANHB;
      // EXTERNAL DGET06, ZLANGE, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZCOPY, ZERRPO, ZGET04, ZLACPY, ZLAIPD, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPBCON, ZPBRFS, ZPBT01, ZPBT02, ZPBT05, ZPBTRF, ZPBTRS, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
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
      // ..
      // .. Executable Statements ..
*
      // Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'PB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      // Test the error exits
*
      IF( TSTERR ) CALL ZERRPO( PATH, NOUT )
      INFOT = 0
      KDVAL( 1 ) = 0
*
      // Do for each value of N in NVAL
*
      DO 90 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
*
         // Set limits on the number of loop iterations.
*
         NKD = MAX( 1, MIN( N, 4 ) )
         NIMAT = NTYPES
         IF( N.EQ.0 ) NIMAT = 1
*
         KDVAL( 2 ) = N + ( N+1 ) / 4
         KDVAL( 3 ) = ( 3*N-1 ) / 4
         KDVAL( 4 ) = ( N+1 ) / 4
*
         DO 80 IKD = 1, NKD
*
            // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            // makes it easier to skip redundant values for small values
            // of N.
*
            KD = KDVAL( IKD )
            LDAB = KD + 1
*
            // Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 70 IUPLO = 1, 2
               KOFF = 1
               IF( IUPLO.EQ.1 ) THEN
                  UPLO = 'U'
                  KOFF = MAX( 1, KD+2-N )
                  PACKIT = 'Q'
               ELSE
                  UPLO = 'L'
                  PACKIT = 'B'
               END IF
*
               DO 60 IMAT = 1, NIMAT
*
                  // Do the tests only if DOTYPE( IMAT ) is true.
*
                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 60
*
                  // Skip types 2, 3, or 4 if the matrix size is too small.
*
                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 ) GO TO 60
*
                  IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) THEN
*
                     // Set up parameters with ZLATB4 and generate a test
                     // matrix with ZLATMS.
*
                     CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
                     SRNAMT = 'ZLATMS'
                     CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD, PACKIT, A( KOFF ), LDAB, WORK, INFO )
*
                     // Check error code from ZLATMS.
*
                     IF( INFO.NE.0 ) THEN
                        CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 60
                     END IF
                  ELSE IF( IZERO.GT.0 ) THEN
*
                     // Use the same matrix for types 3 and 4 as for type
                     // 2 by copying back the zeroed out column,
*
                     IW = 2*LDA + 1
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL ZCOPY( IZERO-I1, WORK( IW ), 1, A( IOFF-IZERO+I1 ), 1 )
                        IW = IW + IZERO - I1
                        CALL ZCOPY( I2-IZERO+1, WORK( IW ), 1, A( IOFF ), MAX( LDAB-1, 1 ) )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL ZCOPY( IZERO-I1, WORK( IW ), 1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ) )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL ZCOPY( I2-IZERO+1, WORK( IW ), 1, A( IOFF ), 1 )
                     END IF
                  END IF
*
                  // For types 2-4, zero one row and column of the matrix
                 t // o test that INFO is returned correctly.
*
                  IZERO = 0
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.2 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.3 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
                     // Save the zeroed out row and column in WORK(*,3)
*
                     IW = 2*LDA
                     DO 20 I = 1, MIN( 2*KD+1, N )
                        WORK( IW+I ) = ZERO
   20                CONTINUE
                     IW = IW + 1
                     I1 = MAX( IZERO-KD, 1 )
                     I2 = MIN( IZERO+KD, N )
*
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL ZSWAP( IZERO-I1, A( IOFF-IZERO+I1 ), 1, WORK( IW ), 1 )
                        IW = IW + IZERO - I1
                        CALL ZSWAP( I2-IZERO+1, A( IOFF ), MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL ZSWAP( IZERO-I1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL ZSWAP( I2-IZERO+1, A( IOFF ), 1, WORK( IW ), 1 )
                     END IF
                  END IF
*
                  // Set the imaginary part of the diagonals.
*
                  IF( IUPLO.EQ.1 ) THEN
                     CALL ZLAIPD( N, A( KD+1 ), LDAB, 0 )
                  ELSE
                     CALL ZLAIPD( N, A( 1 ), LDAB, 0 )
                  END IF
*
                  // Do for each value of NB in NBVAL
*
                  DO 50 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
*
                     // Compute the L*L' or U'*U factorization of the band
                     // matrix.
*
                     CALL ZLACPY( 'Full', KD+1, N, A, LDAB, AFAC, LDAB )
                     SRNAMT = 'ZPBTRF'
                     CALL ZPBTRF( UPLO, N, KD, AFAC, LDAB, INFO )
*
                     // Check error code from ZPBTRF.
*
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'ZPBTRF', INFO, IZERO, UPLO, N, N, KD, KD, NB, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 50
                     END IF
*
                     // Skip the tests if INFO is not 0.
*
                     IF( INFO.NE.0 ) GO TO 50
*
*+    TEST 1
                     // Reconstruct matrix from factors and compute
                     // residual.
*
                     CALL ZLACPY( 'Full', KD+1, N, AFAC, LDAB, AINV, LDAB )                      CALL ZPBT01( UPLO, N, KD, A, LDAB, AINV, LDAB, RWORK, RESULT( 1 ) )
*
                     // Print the test ratio if it is .GE. THRESH.
*
                     IF( RESULT( 1 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, KD, NB, IMAT, 1, RESULT( 1 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
*
                     // Only do other tests if this is the first blocksize.
*
                     IF( INB.GT.1 ) GO TO 50
*
                     // Form the inverse of A so we can get a good estimate
                     // of RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     CALL ZLASET( 'Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), AINV, LDA )
                     SRNAMT = 'ZPBTRS'
                     CALL ZPBTRS( UPLO, N, KD, N, AFAC, LDAB, AINV, LDA, INFO )
*
                     // Compute RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     ANORM = ZLANHB( '1', UPLO, N, KD, A, LDAB, RWORK )
                     AINVNM = ZLANGE( '1', N, N, AINV, LDA, RWORK )
                     IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
*
                     DO 40 IRHS = 1, NNS
                        NRHS = NSVAL( IRHS )
*
*+    TEST 2
                     // Solve and compute residual for A * X = B.
*
                        SRNAMT = 'ZLARHS'
                        CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A, LDAB, XACT, LDA, B, LDA, ISEED, INFO )
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                        SRNAMT = 'ZPBTRS'
                        CALL ZPBTRS( UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA, INFO )
*
                     // Check error code from ZPBTRS.
*
                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZPBTRS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL ZPBT02( UPLO, N, KD, NRHS, A, LDAB, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )
*
*+    TEST 3
                     // Check solution from generated exact solution.
*
                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
*
*+    TESTS 4, 5, and 6
                     // Use iterative refinement to improve the solution.
*
                        SRNAMT = 'ZPBRFS'
                        CALL ZPBRFS( UPLO, N, KD, NRHS, A, LDAB, AFAC, LDAB, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )
*
                     // Check error code from ZPBRFS.
*
                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZPBRFS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )                         CALL ZPBT05( UPLO, N, KD, NRHS, A, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )
*
                        // Print information about the tests that did not
                        // pass the threshold.
*
                        DO 30 K = 2, 6
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, N, KD, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   30                   CONTINUE
                        NRUN = NRUN + 5
   40                CONTINUE
*
*+    TEST 7
                     // Get an estimate of RCOND = 1/CNDNUM.
*
                     SRNAMT = 'ZPBCON'
                     CALL ZPBCON( UPLO, N, KD, AFAC, LDAB, ANORM, RCOND, WORK, RWORK, INFO )
*
                     // Check error code from ZPBCON.
*
                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZPBCON', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT )
*
                     RESULT( 7 ) = DGET06( RCOND, RCONDC )
*
                     // Print the test ratio if it is .GE. THRESH.
*
                     IF( RESULT( 7 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )UPLO, N, KD, IMAT, 7, RESULT( 7 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      // Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NB=', I4,
     $      ', type ', I2, ', test ', I2, ', ratio= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I3,
     $      ', type ', I2, ', test(', I2, ') = ', G12.5 )
 9997 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ',', 10X,
     $      ' type ', I2, ', test(', I2, ') = ', G12.5 )
      RETURN
*
      // End of ZCHKPB
*
      END
