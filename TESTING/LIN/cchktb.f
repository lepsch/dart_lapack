      SUBROUTINE CCHKTB( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, AB, AINV, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      REAL               RWORK( * )
      COMPLEX            AB( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPE1, NTYPES;
      const              NTYPE1 = 9, NTYPES = 17 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IK, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, J, K, KD, LDA, LDAB, N, NERRS, NFAIL, NIMAT, NIMAT2, NK, NRHS, NRUN;
      REAL               AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, SCALE
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANTB, CLANTR
      // EXTERNAL LSAME, CLANTB, CLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRTR, CGET04, CLACPY, CLARHS, CLASET, CLATBS, CLATTB, CTBCON, CTBRFS, CTBSV, CTBT02, CTBT03, CTBT05, CTBT06, CTBTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'TB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL CERRTR( PATH, NOUT )
      INFOT = 0

      DO 140 IN = 1, NN

         // Do for each value of N in NVAL

         N = NVAL( IN )
         LDA = MAX( 1, N )
         XTYPE = 'N'
         NIMAT = NTYPE1
         NIMAT2 = NTYPES
         if ( N.LE.0 ) {
            NIMAT = 1
            NIMAT2 = NTYPE1 + 1
         }

         NK = MIN( N+1, 4 )
         DO 130 IK = 1, NK

            // Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
            // it easier to skip redundant values for small values of N.

            if ( IK.EQ.1 ) {
               KD = 0
            } else if ( IK.EQ.2 ) {
               KD = MAX( N, 0 )
            } else if ( IK.EQ.3 ) {
               KD = ( 3*N-1 ) / 4
            } else if ( IK.EQ.4 ) {
               KD = ( N+1 ) / 4
            }
            LDAB = KD + 1

            DO 90 IMAT = 1, NIMAT

               // Do the tests only if DOTYPE( IMAT ) is true.

               IF( .NOT.DOTYPE( IMAT ) ) GO TO 90

               DO 80 IUPLO = 1, 2

                  // Do first for UPLO = 'U', then for UPLO = 'L'

                  UPLO = UPLOS( IUPLO )

                  // Call CLATTB to generate a triangular test matrix.

                  SRNAMT = 'CLATTB'
                  CALL CLATTB( IMAT, UPLO, 'No transpose', DIAG, ISEED, N, KD, AB, LDAB, X, WORK, RWORK, INFO )

                  // Set IDIAG = 1 for non-unit matrices, 2 for unit.

                  if ( LSAME( DIAG, 'N' ) ) {
                     IDIAG = 1
                  } else {
                     IDIAG = 2
                  }

                  // Form the inverse of A so we can get a good estimate
                  // of RCONDC = 1/(norm(A) * norm(inv(A))).

                  CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), AINV, LDA )
                  if ( LSAME( UPLO, 'U' ) ) {
                     DO 20 J = 1, N
                        CALL CTBSV( UPLO, 'No transpose', DIAG, J, KD, AB, LDAB, AINV( ( J-1 )*LDA+1 ), 1 )
   20                CONTINUE
                  } else {
                     DO 30 J = 1, N
                        CALL CTBSV( UPLO, 'No transpose', DIAG, N-J+1, KD, AB( ( J-1 )*LDAB+1 ), LDAB, AINV( ( J-1 )*LDA+J ), 1 )
   30                CONTINUE
                  }

                  // Compute the 1-norm condition number of A.

                  ANORM = CLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, RWORK )                   AINVNM = CLANTR( '1', UPLO, DIAG, N, N, AINV, LDA, RWORK )
                  if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                     RCONDO = ONE
                  } else {
                     RCONDO = ( ONE / ANORM ) / AINVNM
                  }

                  // Compute the infinity-norm condition number of A.

                  ANORM = CLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, RWORK )                   AINVNM = CLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK )
                  if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                     RCONDI = ONE
                  } else {
                     RCONDI = ( ONE / ANORM ) / AINVNM
                  }

                  DO 60 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'

                     DO 50 ITRAN = 1, NTRAN

                     // Do for op(A) = A, A**T, or A**H.

                        TRANS = TRANSS( ITRAN )
                        if ( ITRAN.EQ.1 ) {
                           NORM = 'O'
                           RCONDC = RCONDO
                        } else {
                           NORM = 'I'
                           RCONDC = RCONDI
                        }

*+    TEST 1
                     // Solve and compute residual for op(A)*x = b.

                        SRNAMT = 'CLARHS'
                        CALL CLARHS( PATH, XTYPE, UPLO, TRANS, N, N, KD, IDIAG, NRHS, AB, LDAB, XACT, LDA, B, LDA, ISEED, INFO )
                        XTYPE = 'C'
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                        SRNAMT = 'CTBTRS'
                        CALL CTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDA, INFO )

                     // Check error code from CTBTRS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'CTBTRS', INFO, 0, UPLO // TRANS // DIAG, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        CALL CTBT02( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDA, B, LDA, WORK, RWORK, RESULT( 1 ) )

*+    TEST 2
                     // Check solution from generated exact solution.

                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 2 ) )

*+    TESTS 3, 4, and 5
                     // Use iterative refinement to improve the solution
                     // and compute error bounds.

                        SRNAMT = 'CTBRFS'
                        CALL CTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                     // Check error code from CTBRFS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'CTBRFS', INFO, 0, UPLO // TRANS // DIAG, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )                         CALL CTBT05( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 40 K = 1, 5
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )UPLO, TRANS, DIAG, N, KD, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
   40                   CONTINUE
                        NRUN = NRUN + 5
   50                CONTINUE
   60             CONTINUE

*+    TEST 6
                     // Get an estimate of RCOND = 1/CNDNUM.

                  DO 70 ITRAN = 1, 2
                     if ( ITRAN.EQ.1 ) {
                        NORM = 'O'
                        RCONDC = RCONDO
                     } else {
                        NORM = 'I'
                        RCONDC = RCONDI
                     }
                     SRNAMT = 'CTBCON'
                     CALL CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO )

                     // Check error code from CTBCON.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'CTBCON', INFO, 0, NORM // UPLO // DIAG, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL CTBT06( RCOND, RCONDC, UPLO, DIAG, N, KD, AB, LDAB, RWORK, RESULT( 6 ) )

                     // Print the test ratio if it is .GE. THRESH.

                     if ( RESULT( 6 ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 ) 'CTBCON', NORM, UPLO, DIAG, N, KD, IMAT, 6, RESULT( 6 )
                        NFAIL = NFAIL + 1
                     }
                     NRUN = NRUN + 1
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE

            // Use pathological test matrices to test CLATBS.

            DO 120 IMAT = NTYPE1 + 1, NIMAT2

               // Do the tests only if DOTYPE( IMAT ) is true.

               IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

               DO 110 IUPLO = 1, 2

                  // Do first for UPLO = 'U', then for UPLO = 'L'

                  UPLO = UPLOS( IUPLO )
                  DO 100 ITRAN = 1, NTRAN

                     // Do for op(A) = A, A**T, and A**H.

                     TRANS = TRANSS( ITRAN )

                     // Call CLATTB to generate a triangular test matrix.

                     SRNAMT = 'CLATTB'
                     CALL CLATTB( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB, LDAB, X, WORK, RWORK, INFO )

*+    TEST 7
                     // Solve the system op(A)*x = b

                     SRNAMT = 'CLATBS'
                     CALL CCOPY( N, X, 1, B, 1 )
                     CALL CLATBS( UPLO, TRANS, DIAG, 'N', N, KD, AB, LDAB, B, SCALE, RWORK, INFO )

                     // Check error code from CLATBS.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'CLATBS', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL CTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 7 ) )

*+    TEST 8
                     // Solve op(A)*x = b again with NORMIN = 'Y'.

                     CALL CCOPY( N, X, 1, B, 1 )
                     CALL CLATBS( UPLO, TRANS, DIAG, 'Y', N, KD, AB, LDAB, B, SCALE, RWORK, INFO )

                     // Check error code from CLATBS.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'CLATBS', INFO, 0, UPLO // TRANS // DIAG // 'Y', N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL CTBT03( UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( RESULT( 7 ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )'CLATBS', UPLO, TRANS, DIAG, 'N', N, KD, IMAT, 7, RESULT( 7 )
                        NFAIL = NFAIL + 1
                     }
                     if ( RESULT( 8 ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )'CLATBS', UPLO, TRANS, DIAG, 'Y', N, KD, IMAT, 8, RESULT( 8 )
                        NFAIL = NFAIL + 1
                     }
                     NRUN = NRUN + 2
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I5, ', type ', I2, ', test(', I2, ')=', G12.5 )
 9998 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''',', I5, ',', I5, ',  ... ), type ', I2, ', test(', I2, ')=', G12.5 )
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''', A1, ''',', I5, ',', I5, ', ...  ),  type ', I2, ', test(', I1, ')=', G12.5 )
      RETURN

      // End of CCHKTB

      }
