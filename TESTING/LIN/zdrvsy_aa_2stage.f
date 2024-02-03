      SUBROUTINE ZDRVSY_AA_2STAGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 3 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, TYPE, UPLO, XTYPE;
      String             MATPATH, PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             ANORM, CNDNUM;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANSY;
      // EXTERNAL DGET06, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGET04, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZSYSV_AA_2STAGE, ZSYT01_AA, ZSYT02, ZSYTRF_AA_2STAGE
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
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      // Test path

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'S2'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Zomplex precision'
      MATPATH( 2: 3 ) = 'SY'

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

      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

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

               // Set up parameters with ZLATB4 for the matrix generator
               // based on the type of matrix to be generated.

              CALL ZLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               // Generate a matrix with ZLATMS.

                  SRNAMT = 'ZLATMS'
                  CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )

                  // Check error code from ZLATMS and handle error.

                  if ( INFO.NE.0 ) {
                     CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 160
                  }

                  // For types 3-6, zero one or more rows and columns of
                  // the matrix to test that INFO is returned correctly.

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
                              A( IOFF+I ) = CZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        } else {
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = CZERO
   50                      CONTINUE
                        }
                     } else {
                        IOFF = 0
                        if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = CZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                           IZERO = 1
                        } else {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = CZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        }
                     }
                  } else {
                     IZERO = 0
                  }

                  // End generate the test matrix A.


               DO 150 IFACT = 1, NFACT

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                  XTYPE = 'C'

                  // --- Test ZSYSV_AA_2STAGE  ---

                  if ( IFACT.EQ.2 ) {
                     CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     // Factor the matrix and solve the system using ZSYSV_AA.

                     SRNAMT = 'ZSYSV_AA_2STAGE '
                     LWORK = MIN(N*NB, 3*NMAX*NMAX)
                     CALL ZSYSV_AA_2STAGE( UPLO, N, NRHS, AFAC, LDA, AINV, (3*NB+1)*N, IWORK, IWORK( 1+N ), X, LDA, WORK, LWORK, INFO )

                     // Adjust the expected value of INFO to account for
                     // pivoting.

                     if ( IZERO.GT.0 ) {
                        J = 1
                        K = IZERO
  100                   CONTINUE
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

                     // Check error code from ZSYSV_AA_2STAGE .

                     if ( INFO.NE.K ) {
                        CALL ALAERH( PATH, 'ZSYSV_AA_2STAGE', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

                     // Compute residual of the computed solution.

                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) )

                     // Reconstruct matrix from factors and compute
                     // residual.

                      // CALL ZSY01_AA( UPLO, N, A, LDA, AFAC, LDA,
      // $                                  IWORK, AINV, LDA, RWORK,
      // $                                  RESULT( 2 ) )
                      // NT = 2
                     NT = 1

                     // Print information about the tests that did not pass
                     // the threshold.

                     DO 110 K = 1, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'ZSYSV_AA_2STAGE ', UPLO, N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        }
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  }

  150          CONTINUE

  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of ZDRVSY_AA_2STAGE

      }
