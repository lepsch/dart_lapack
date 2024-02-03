      SUBROUTINE SCHKPS( DOTYPE, NN, NVAL, NNB, NBVAL, NRANK, RANKVAL, THRESH, TSTERR, NMAX, A, AFAC, PERM, PIV, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               THRESH
      int                NMAX, NN, NNB, NOUT, NRANK;
      bool               TSTERR;
      // ..
      // .. Array Arguments ..
      REAL               A( * ), AFAC( * ), PERM( * ), RWORK( * ), WORK( * )
      int                NBVAL( * ), NVAL( * ), PIV( * ), RANKVAL( * );
      bool               DOTYPE( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      int                NTYPES;
      PARAMETER          ( NTYPES = 9 )
      // ..
      // .. Local Scalars ..
      REAL               ANORM, CNDNUM, RESULT, TOL
      int                COMPRANK, I, IMAT, IN, INB, INFO, IRANK, IUPLO, IZERO, KL, KU, LDA, MODE, N, NB, NERRS, NFAIL, NIMAT, NRUN, RANK, RANKDIFF;
      String             DIST, TYPE, UPLO;
      String             PATH;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      String             UPLOS( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRPS, SLACPY, SLATB5, SLATMT, SPST01, SPSTRF, XLAENV
      // ..
      // .. Scalars in Common ..
      int                INFOT, NUNIT;
      bool               LERR, OK;
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, CEILING
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Single Precision'
      PATH( 2: 3 ) = 'PS'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 100 I = 1, 4
         ISEED( I ) = ISEEDY( I )
  100 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL SERRPS( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )

      // Do for each value of N in NVAL

      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         IZERO = 0
         DO 140 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 140

               // Do for each value of RANK in RANKVAL

            DO 130 IRANK = 1, NRANK

               // Only repeat test 3 to 5 for different ranks
               // Other tests use full rank

               IF( ( IMAT.LT.3 .OR. IMAT.GT.5 ) .AND. IRANK.GT.1 ) GO TO 130

               RANK = CEILING( ( N * REAL( RANKVAL( IRANK ) ) ) / 100.E+0 )


            // Do first for UPLO = 'U', then for UPLO = 'L'

               DO 120 IUPLO = 1, 2
                  UPLO = UPLOS( IUPLO )

               // Set up parameters with SLATB5 and generate a test matrix
               // with SLATMT.

                  CALL SLATB5( PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

                  SRNAMT = 'SLATMT'
                  CALL SLATMT( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, RANK, KL, KU, UPLO, A, LDA, WORK, INFO )

               // Check error code from SLATMT.

                  IF( INFO.NE.0 ) THEN
                    CALL ALAERH( PATH, 'SLATMT', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 120
                  END IF

               // Do for each value of NB in NBVAL

                  DO 110 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )

                  // Compute the pivoted L*L' or U'*U factorization
                  // of the matrix.

                     CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     SRNAMT = 'SPSTRF'

                  // Use default tolerance

                     TOL = -ONE
                     CALL SPSTRF( UPLO, N, AFAC, LDA, PIV, COMPRANK, TOL, WORK, INFO )

                  // Check error code from SPSTRF.

                     IF( (INFO.LT.IZERO) .OR.(INFO.NE.IZERO.AND.RANK.EQ.N) .OR.(INFO.LE.IZERO.AND.RANK.LT.N) ) THEN                         CALL ALAERH( PATH, 'SPSTRF', INFO, IZERO, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 110
                     END IF

                  // Skip the test if INFO is not 0.

                     IF( INFO.NE.0 ) GO TO 110

                  // Reconstruct matrix from factors and compute residual.

                  // PERM holds permuted L*L^T or U^T*U

                     CALL SPST01( UPLO, N, A, LDA, AFAC, LDA, PERM, LDA, PIV, RWORK, RESULT, COMPRANK )

                  // Print information about the tests that did not pass
                 t // he threshold or where computed rank was not RANK.

                     IF( N.EQ.0 ) COMPRANK = 0
                     RANKDIFF = RANK - COMPRANK
                     IF( RESULT.GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, RANK, RANKDIFF, NB, IMAT, RESULT
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
  110             CONTINUE

  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', RANK =', I3,
     $      ', Diff =', I5, ', NB =', I4, ', type ', I2, ', Ratio =',
     $      G12.5 )
      RETURN

      // End of SCHKPS

      END
