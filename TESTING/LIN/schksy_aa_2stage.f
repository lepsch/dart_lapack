      SUBROUTINE SCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      bool         TSTERR;
      int          NN, NNB, NNS, NMAX, NOUT;
      REAL         THRESH
      // ..
      // .. Array Arguments ..
      bool         DOTYPE( * );
      int          IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      REAL         A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL         ZERO
      const        ZERO = 0.0E+0 ;
      int          NTYPES;
      const        NTYPES = 10 ;
      int          NTESTS;
      const        NTESTS = 9 ;
      // ..
      // .. Local Scalars ..
      bool         ZEROT;
      String       DIST, TYPE, UPLO, XTYPE;
      String       PATH, MATPATH;
      int          I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      REAL         ANORM, CNDNUM
      // ..
      // .. Local Arrays ..
      String       UPLOS( 2 );
      int          ISEED( 4 ), ISEEDY( 4 );
      REAL         RESULT( NTESTS )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRSY, SLACPY, SLARHS, SLATB4, SLATMS, SPOT02, SSYT01_AA, SSYTRF_AA_2STAGE, SSYTRS_AA_2STAGE, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool         LERR, OK;
      String       SRNAMT;
      int          INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON       / INFOC / INFOT, NUNIT, OK, LERR
      COMMON       / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA         ISEEDY / 1988, 1989, 1990, 1991 /
      DATA         UPLOS / 'U', 'L' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.


      // Test path

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'S2'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Single precision'
      MATPATH( 2: 3 ) = 'SY'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL SERRSY( PATH, NOUT )
      INFOT = 0

      // Set the minimum block size for which the block routine should
      // be used, which will be later returned by ILAENV

      CALL XLAENV( 2, 2 )

      // Do for each value of N in NVAL

      DO 180 IN = 1, NN
         N = NVAL( IN )
         if ( N .GT. NMAX ) {
            NFAIL = NFAIL + 1
            WRITE(NOUT, 9995) 'M ', N, NMAX
            GO TO 180
         }
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         IZERO = 0

         // Do for each value of matrix type IMAT

         DO 170 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 170

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 160 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )

               // Begin generate the test matrix A.


               // Set up parameters with SLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               CALL SLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               // Generate a matrix with SLATMS.

               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )

               // Check error code from SLATMS and handle error.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                     // Skip all tests for this generated matrix

                  GO TO 160
               }

               // For matrix types 3-6, zero one or more rows and
               // columns of the matrix to test that INFO is returned
               // correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }

                  if ( IMAT.LT.6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO.EQ.1 ) {
                        IOFF = ( IZERO-1 )*LDA
                        DO 20 I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                        IOFF = IOFF + IZERO
                        DO 30 I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   30                   CONTINUE
                     } else {
                        IOFF = IZERO
                        DO 40 I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   40                   CONTINUE
                        IOFF = IOFF - IZERO
                        DO 50 I = IZERO, N
                           A( IOFF+I ) = ZERO
   50                   CONTINUE
                     }
                  } else {
                     if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                        IOFF = 0
                        DO 70 J = 1, N
                           I2 = MIN( J, IZERO )
                           DO 60 I = 1, I2
                              A( IOFF+I ) = ZERO
   60                      CONTINUE
                           IOFF = IOFF + LDA
   70                   CONTINUE
                        IZERO = 1
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        IOFF = 0
                        DO 90 J = 1, N
                           I1 = MAX( J, IZERO )
                           DO 80 I = I1, N
                              A( IOFF+I ) = ZERO
   80                      CONTINUE
                           IOFF = IOFF + LDA
   90                   CONTINUE
                     }
                  }
               } else {
                  IZERO = 0
               }

               // End generate the test matrix A.

               // Do for each value of NB in NBVAL

               DO 150 INB = 1, NNB

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                 t // he block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                  SRNAMT = 'SSYTRF_AA_2STAGE'
                  LWORK = MIN( MAX( 1, N*NB ), 3*NMAX*NMAX )
                  CALL SSYTRF_AA_2STAGE( UPLO, N, AFAC, LDA,  AINV, MAX( 1, (3*NB+1)*N ), IWORK, IWORK( 1+N ), WORK, LWORK, INFO )

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  if ( IZERO.GT.0 ) {
                     J = 1
                     K = IZERO
  100                CONTINUE
                     if ( J.EQ.K ) {
                        K = IWORK( J )
                     } else if ( IWORK( J ).EQ.K ) {
                        K = J
                     }
                     if ( J.LT.K ) {
                        J = J + 1
                        GO TO 100
                     }
                  } else {
                     K = 0
                  }

                  // Check error code from SSYTRF and handle error.

                  if ( INFO.NE.K ) {
                     CALL ALAERH( PATH, 'SSYTRF_AA_2STAGE', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
                  }

*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                   // CALL SSYT01_AA( UPLO, N, A, LDA, AFAC, LDA, IWORK,
      // $                            AINV, LDA, RWORK, RESULT( 1 ) )
                   // NT = 1
                   NT = 0


                  // Print information about the tests that did not pass
                 t // he threshold.

                  DO 110 K = 1, NT
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
  110             CONTINUE
                  NRUN = NRUN + NT

                  // Skip solver test if INFO is not 0.

                  if ( INFO.NE.0 ) {
                     GO TO 140
                  }

                  // Do for each value of NRHS in NSVAL.

                  DO 130 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )

*+    TEST 2 (Using TRS)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                     SRNAMT = 'SLARHS'
                     CALL SLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     SRNAMT = 'SSYTRS_AA_2STAGE'
                     CALL SSYTRS_AA_2STAGE( UPLO, N, NRHS, AFAC, LDA, AINV, (3*NB+1)*N, IWORK, IWORK( 1+N ), X, LDA, INFO )

                     // Check error code from SSYTRS and handle error.

                     if ( INFO.NE.0 ) {
                        if ( IZERO.EQ.0 ) {
                           CALL ALAERH( PATH, 'SSYTRS_AA_2STAGE', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        }
                     } else {
                        CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )

                        // Compute the residual for the solution

                        CALL SPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )


                     // Print information about the tests that did not pass
                    t // he threshold.

                        DO 120 K = 2, 2
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
  120                   CONTINUE
                     }
                     NRUN = NRUN + 1

                  // End do for each value of NRHS in NSVAL.

  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=',
     $      I6 )
      RETURN

      // End of SCHKSY_AA_2STAGE

      }
