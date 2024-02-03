      SUBROUTINE SDRVSY_AA_2STAGE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      REAL               RWORK( * )
      REAL               A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
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
      REAL               ANORM, CNDNUM
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      REAL               SLANSY, SGET06
      // EXTERNAL SLANSY, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, SERRVX, SLACPY, SLARHS, SLATB4, SLATMS, SSYSV_AA_2STAGE, SSYT01_AA, SPOT02, SSYTRF_AA_2STAGE
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

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'S2'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Single precision'
      MATPATH( 2: 3 ) = 'SY'

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL SERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 170

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO )

               // Begin generate the test matrix A.

               // Set up parameters with SLATB4 for the matrix generator
               // based on the type of matrix to be generated.

              slatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with SLATMS.

                  SRNAMT = 'SLATMS'
                  slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

                  // Check error code from SLATMS and handle error.

                  if ( INFO.NE.0 ) {
                     alaerh(PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
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
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           for (I = IZERO; I <= N; I++) { // 30
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        } else {
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           for (I = IZERO; I <= N; I++) { // 50
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        }
                     } else {
                        IOFF = 0
                        if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                           for (J = 1; J <= N; J++) { // 70
                              I2 = MIN( J, IZERO )
                              for (I = 1; I <= I2; I++) { // 60
                                 A( IOFF+I ) = ZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                           IZERO = 1
                        } else {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           for (J = 1; J <= N; J++) { // 90
                              I1 = MAX( J, IZERO )
                              for (I = I1; I <= N; I++) { // 80
                                 A( IOFF+I ) = ZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        }
                     }
                  } else {
                     IZERO = 0
                  }

                  // End generate the test matrix A.


               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'SLARHS'
                  slarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C'

                  // --- Test SSYSV_AA_2STAGE  ---

                  if ( IFACT.EQ.2 ) {
                     slacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     slacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using SSYSV_AA.

                     SRNAMT = 'SSYSV_AA_2STAGE '
                     LWORK = MIN( MAX( 1, N*NB ), 3*NMAX*NMAX )
                     ssysv_aa_2stage(UPLO, N, NRHS, AFAC, LDA, AINV, MAX( 1, (3*NB+1)*N ), IWORK, IWORK( 1+N ), X, LDA, WORK, LWORK, INFO );

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

                     // Check error code from SSYSV_AA .

                     if ( INFO.NE.K ) {
                        alaerh(PATH, 'SSYSV_AA', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

                     // Compute residual of the computed solution.

                     slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     spot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) );

                     // Reconstruct matrix from factors and compute
                     // residual.

                      // CALL SSY01_AA( UPLO, N, A, LDA, AFAC, LDA,
      // $                                  IWORK, AINV, LDA, RWORK,
      // $                                  RESULT( 2 ) )
                      // NT = 2
                     NT = 1

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'SSYSV_AA ', UPLO, N, IMAT, K, RESULT( K )
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

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of SDRVSY_AA_2STAGE

      }
