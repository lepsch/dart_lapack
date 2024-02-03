      SUBROUTINE SCHKGB( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, A, LA, AFAC, LAFAC, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                LA, LAFAC, NM, NN, NNB, NNS, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      REAL               A( * ), AFAC( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 7 ;
      int                NBW, NTRAN;
      const              NBW = 4, NTRAN = 3 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, I1, I2, IKL, IKU, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, J, K, KL, KOFF, KU, LDA, LDAFAC, LDB, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NKL, NKU, NRHS, NRUN;
      REAL               AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, RCOND, RCONDC, RCONDI, RCONDO
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 ), KLVAL( NBW ), KUVAL( NBW );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      REAL               SGET06, SLANGB, SLANGE
      // EXTERNAL SGET06, SLANGB, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SCOPY, SERRGE, SGBCON, SGBRFS, SGBT01, SGBT02, SGBT05, SGBTRF, SGBTRS, SGET04, SLACPY, SLARHS, SLASET, SLATB4, SLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
      DATA               ISEEDY / 1988, 1989, 1990, 1991 / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'GB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL SERRGE( PATH, NOUT )
      INFOT = 0
      xlaenv(2, 2 );

      // Initialize the first value for the lower and upper bandwidths.

      KLVAL( 1 ) = 0
      KUVAL( 1 ) = 0

      // Do for each value of M in MVAL

      DO 160 IM = 1, NM
         M = MVAL( IM )

         // Set values to use for the lower bandwidth.

         KLVAL( 2 ) = M + ( M+1 ) / 4

         // KLVAL( 2 ) = MAX( M-1, 0 )

         KLVAL( 3 ) = ( 3*M-1 ) / 4
         KLVAL( 4 ) = ( M+1 ) / 4

         // Do for each value of N in NVAL

         DO 150 IN = 1, NN
            N = NVAL( IN )
            XTYPE = 'N'

            // Set values to use for the upper bandwidth.

            KUVAL( 2 ) = N + ( N+1 ) / 4

            // KUVAL( 2 ) = MAX( N-1, 0 )

            KUVAL( 3 ) = ( 3*N-1 ) / 4
            KUVAL( 4 ) = ( N+1 ) / 4

            // Set limits on the number of loop iterations.

            NKL = MIN( M+1, 4 )
            IF( N.EQ.0 ) NKL = 2
            NKU = MIN( N+1, 4 )
            IF( M.EQ.0 ) NKU = 2
            NIMAT = NTYPES
            IF( M.LE.0 .OR. N.LE.0 ) NIMAT = 1

            DO 140 IKL = 1, NKL

               // Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
               // order makes it easier to skip redundant values for small
               // values of M.

               KL = KLVAL( IKL )
               DO 130 IKU = 1, NKU

                  // Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
                  // order makes it easier to skip redundant values for
                  // small values of N.

                  KU = KUVAL( IKU )

                  // Check that A and AFAC are big enough to generate this
                  // matrix.

                  LDA = KL + KU + 1
                  LDAFAC = 2*KL + KU + 1
                  if ( ( LDA*N ).GT.LA .OR. ( LDAFAC*N ).GT.LAFAC ) {
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
                     if ( N*( KL+KU+1 ).GT.LA ) {
                        WRITE( NOUT, FMT = 9999 )LA, M, N, KL, KU, N*( KL+KU+1 )
                        NERRS = NERRS + 1
                     }
                     if ( N*( 2*KL+KU+1 ).GT.LAFAC ) {
                        WRITE( NOUT, FMT = 9998 )LAFAC, M, N, KL, KU, N*( 2*KL+KU+1 )
                        NERRS = NERRS + 1
                     }
                     GO TO 130
                  }

                  DO 120 IMAT = 1, NIMAT

                     // Do the tests only if DOTYPE( IMAT ) is true.

                     IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

                     // Skip types 2, 3, or 4 if the matrix size is too
                     // small.

                     ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                     IF( ZEROT .AND. N.LT.IMAT-1 ) GO TO 120

                     if ( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) {

                        // Set up parameters with SLATB4 and generate a
                        // test matrix with SLATMS.

                        slatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                        KOFF = MAX( 1, KU+2-N )
                        DO 20 I = 1, KOFF - 1
                           A( I ) = ZERO
   20                   CONTINUE
                        SRNAMT = 'SLATMS'
                        slatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A( KOFF ), LDA, WORK, INFO );

                        // Check the error code from SLATMS.

                        if ( INFO.NE.0 ) {
                           alaerh(PATH, 'SLATMS', INFO, 0, ' ', M, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 120
                        }
                     } else if ( IZERO.GT.0 ) {

                        // Use the same matrix for types 3 and 4 as for
                        // type 2 by copying back the zeroed out column.

                        scopy(I2-I1+1, B, 1, A( IOFF+I1 ), 1 );
                     }

                     // For types 2, 3, and 4, zero one or more columns of
                     // the matrix to test that INFO is returned correctly.

                     IZERO = 0
                     if ( ZEROT ) {
                        if ( IMAT.EQ.2 ) {
                           IZERO = 1
                        } else if ( IMAT.EQ.3 ) {
                           IZERO = MIN( M, N )
                        } else {
                           IZERO = MIN( M, N ) / 2 + 1
                        }
                        IOFF = ( IZERO-1 )*LDA
                        if ( IMAT.LT.4 ) {

                           // Store the column to be zeroed out in B.

                           I1 = MAX( 1, KU+2-IZERO )
                           I2 = MIN( KL+KU+1, KU+1+( M-IZERO ) )
                           scopy(I2-I1+1, A( IOFF+I1 ), 1, B, 1 );

                           DO 30 I = I1, I2
                              A( IOFF+I ) = ZERO
   30                      CONTINUE
                        } else {
                           DO 50 J = IZERO, N
                              DO 40 I = MAX( 1, KU+2-J ), MIN( KL+KU+1, KU+1+( M-J ) )
                                 A( IOFF+I ) = ZERO
   40                         CONTINUE
                              IOFF = IOFF + LDA
   50                      CONTINUE
                        }
                     }

                     // These lines, if used in place of the calls in the
                     // loop over INB, cause the code to bomb on a Sun
                     // SPARCstation.

                      // ANORMO = SLANGB( 'O', N, KL, KU, A, LDA, RWORK )
                      // ANORMI = SLANGB( 'I', N, KL, KU, A, LDA, RWORK )

                     // Do for each blocksize in NBVAL

                     DO 110 INB = 1, NNB
                        NB = NBVAL( INB )
                        xlaenv(1, NB );

                        // Compute the LU factorization of the band matrix.

                        IF( M.GT.0 .AND. N.GT.0 ) CALL SLACPY( 'Full', KL+KU+1, N, A, LDA, AFAC( KL+1 ), LDAFAC )
                        SRNAMT = 'SGBTRF'
                        sgbtrf(M, N, KL, KU, AFAC, LDAFAC, IWORK, INFO );

                        // Check error code from SGBTRF.

                        IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'SGBTRF', INFO, IZERO, ' ', M, N, KL, KU, NB, IMAT, NFAIL, NERRS, NOUT )
                        TRFCON = .FALSE.

*+    TEST 1
                        // Reconstruct matrix from factors and compute
                        // residual.

                        sgbt01(M, N, KL, KU, A, LDA, AFAC, LDAFAC, IWORK, WORK, RESULT( 1 ) );

                        // Print information about the tests so far that
                        // did not pass the threshold.

                        if ( RESULT( 1 ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9997 )M, N, KL, KU, NB, IMAT, 1, RESULT( 1 )
                           NFAIL = NFAIL + 1
                        }
                        NRUN = NRUN + 1

                        // Skip the remaining tests if this is not the
                        // first block size or if M .ne. N.

                        IF( INB.GT.1 .OR. M.NE.N ) GO TO 110

                        ANORMO = SLANGB( 'O', N, KL, KU, A, LDA, RWORK )
                        ANORMI = SLANGB( 'I', N, KL, KU, A, LDA, RWORK )

                        if ( INFO.EQ.0 ) {

                           // Form the inverse of A so we can get a good
                           // estimate of CNDNUM = norm(A) * norm(inv(A)).

                           LDB = MAX( 1, N )
                           slaset('Full', N, N, ZERO, ONE, WORK, LDB );
                           SRNAMT = 'SGBTRS'
                           sgbtrs('No transpose', N, KL, KU, N, AFAC, LDAFAC, IWORK, WORK, LDB, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = SLANGE( 'O', N, N, WORK, LDB, RWORK )
                           if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDO = ONE
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM
                           }

                           // Compute the infinity-norm condition number of
                           // A.

                           AINVNM = SLANGE( 'I', N, N, WORK, LDB, RWORK )
                           if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDI = ONE
                           } else {
                              RCONDI = ( ONE / ANORMI ) / AINVNM
                           }
                        } else {

                           // Do only the condition estimate if INFO.NE.0.

                           TRFCON = .TRUE.
                           RCONDO = ZERO
                           RCONDI = ZERO
                        }

                        // Skip the solve tests if the matrix is singular.

                        IF( TRFCON ) GO TO 90

                        DO 80 IRHS = 1, NNS
                           NRHS = NSVAL( IRHS )
                           XTYPE = 'N'

                           DO 70 ITRAN = 1, NTRAN
                              TRANS = TRANSS( ITRAN )
                              if ( ITRAN.EQ.1 ) {
                                 RCONDC = RCONDO
                                 NORM = 'O'
                              } else {
                                 RCONDC = RCONDI
                                 NORM = 'I'
                              }

*+    TEST 2:
                              // Solve and compute residual for op(A) * X = B.

                              SRNAMT = 'SLARHS'
                              slarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO );
                              XTYPE = 'C'
                              slacpy('Full', N, NRHS, B, LDB, X, LDB );

                              SRNAMT = 'SGBTRS'
                              sgbtrs(TRANS, N, KL, KU, NRHS, AFAC, LDAFAC, IWORK, X, LDB, INFO );

                              // Check error code from SGBTRS.

                              IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGBTRS', INFO, 0, TRANS, N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )

                              slacpy('Full', N, NRHS, B, LDB, WORK, LDB )                               CALL SGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

*+    TEST 3:
                              // Check solution from generated exact
                              // solution.

                              sget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, 6:
                              // Use iterative refinement to improve the
                              // solution.

                              SRNAMT = 'SGBRFS'
                              sgbrfs(TRANS, N, KL, KU, NRHS, A, LDA, AFAC, LDAFAC, IWORK, B, LDB, X, LDB, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                              // Check error code from SGBRFS.

                              IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGBRFS', INFO, 0, TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT )

                              sget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 4 ) )                               CALL SGBT05( TRANS, N, KL, KU, NRHS, A, LDA, B, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );
                              DO 60 K = 2, 6
                                 if ( RESULT( K ).GE.THRESH ) {
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9996 )TRANS, N, KL, KU, NRHS, IMAT, K, RESULT( K )
                                    NFAIL = NFAIL + 1
                                 }
   60                         CONTINUE
                              NRUN = NRUN + 5
   70                      CONTINUE
   80                   CONTINUE

*+    TEST 7:
                           // Get an estimate of RCOND = 1/CNDNUM.

   90                   CONTINUE
                        DO 100 ITRAN = 1, 2
                           if ( ITRAN.EQ.1 ) {
                              ANORM = ANORMO
                              RCONDC = RCONDO
                              NORM = 'O'
                           } else {
                              ANORM = ANORMI
                              RCONDC = RCONDI
                              NORM = 'I'
                           }
                           SRNAMT = 'SGBCON'
                           sgbcon(NORM, N, KL, KU, AFAC, LDAFAC, IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO );

                              // Check error code from SGBCON.

                           IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGBCON', INFO, 0, NORM, N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )

                           RESULT( 7 ) = SGET06( RCOND, RCONDC )

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( RESULT( 7 ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9995 )NORM, N, KL, KU, IMAT, 7, RESULT( 7 )
                              NFAIL = NFAIL + 1
                           }
                           NRUN = NRUN + 1
  100                   CONTINUE

  110                CONTINUE
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' *** In SCHKGB, LA=', I5, ' is too small for M=', I5, ', N=', I5, ', KL=', I4, ', KU=', I4, / ' ==> Increase LA to at least ', I5 )
 9998 FORMAT( ' *** In SCHKGB, LAFAC=', I5, ' is too small for M=', I5, ', N=', I5, ', KL=', I4, ', KU=', I4, / ' ==> Increase LAFAC to at least ', I5 )
 9997 FORMAT( ' M =', I5, ', N =', I5, ', KL=', I5, ', KU=', I5, ', NB =', I4, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( ' TRANS=''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, ', NRHS=', I3, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( ' NORM =''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, ',', 10X, ' type ', I1, ', test(', I1, ')=', G12.5 )

      RETURN

      // End of SCHKGB

      }
