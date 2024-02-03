      SUBROUTINE DCHKPB( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 7 ;
      int                NBW;
      const              NBW = 4 ;
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
      double             DGET06, DLANGE, DLANSB;
      // EXTERNAL DGET06, DLANGE, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRPO, DGET04, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, DPBCON, DPBRFS, DPBT01, DPBT02, DPBT05, DPBTRF, DPBTRS, DSWAP, XLAENV
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
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /;
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'PB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) CALL DERRPO( PATH, NOUT );
      INFOT = 0;
      xlaenv(2, 2 );
      KDVAL( 1 ) = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 90
         N = NVAL( IN );
         LDA = MAX( N, 1 );
         XTYPE = 'N';

         // Set limits on the number of loop iterations.

         NKD = MAX( 1, MIN( N, 4 ) );
         NIMAT = NTYPES;
         if (N == 0) NIMAT = 1;

         KDVAL( 2 ) = N + ( N+1 ) / 4;
         KDVAL( 3 ) = ( 3*N-1 ) / 4;
         KDVAL( 4 ) = ( N+1 ) / 4;

         for (IKD = 1; IKD <= NKD; IKD++) { // 80

            // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            // makes it easier to skip redundant values for small values
            // of N.

            KD = KDVAL( IKD );
            LDAB = KD + 1;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 70
               KOFF = 1;
               if ( IUPLO == 1 ) {
                  UPLO = 'U';
                  KOFF = MAX( 1, KD+2-N );
                  PACKIT = 'Q';
               } else {
                  UPLO = 'L';
                  PACKIT = 'B';
               }

               for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 60

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  IF( !DOTYPE( IMAT ) ) GO TO 60;

                  // Skip types 2, 3, or 4 if the matrix size is too small.

                  ZEROT = IMAT >= 2 && IMAT <= 4;
                  if (ZEROT && N < IMAT-1) GO TO 60;

                  if ( !ZEROT || !DOTYPE( 1 ) ) {

                     // Set up parameters with DLATB4 and generate a test
                     // matrix with DLATMS.

                     dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                     SRNAMT = 'DLATMS';
                     dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD, PACKIT, A( KOFF ), LDAB, WORK, INFO );

                     // Check error code from DLATMS.

                     if ( INFO != 0 ) {
                        alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 60;
                     }
                  } else if ( IZERO > 0 ) {

                     // Use the same matrix for types 3 and 4 as for type
                     // 2 by copying back the zeroed out column,

                     IW = 2*LDA + 1;
                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        dcopy(IZERO-I1, WORK( IW ), 1, A( IOFF-IZERO+I1 ), 1 );
                        IW = IW + IZERO - I1;
                        dcopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), MAX( LDAB-1, 1 ) );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        dcopy(IZERO-I1, WORK( IW ), 1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ) );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        dcopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), 1 );
                     }
                  }

                  // For types 2-4, zero one row and column of the matrix
                  // to test that INFO is returned correctly.

                  IZERO = 0;
                  if ( ZEROT ) {
                     if ( IMAT == 2 ) {
                        IZERO = 1;
                     } else if ( IMAT == 3 ) {
                        IZERO = N;
                     } else {
                        IZERO = N / 2 + 1;
                     }

                     // Save the zeroed out row and column in WORK(*,3)

                     IW = 2*LDA;
                     DO 20 I = 1, MIN( 2*KD+1, N );
                        WORK( IW+I ) = ZERO;
                     } // 20
                     IW = IW + 1;
                     I1 = MAX( IZERO-KD, 1 );
                     I2 = MIN( IZERO+KD, N );

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        dswap(IZERO-I1, A( IOFF-IZERO+I1 ), 1, WORK( IW ), 1 );
                        IW = IW + IZERO - I1;
                        dswap(I2-IZERO+1, A( IOFF ), MAX( LDAB-1, 1 ), WORK( IW ), 1 );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        dswap(IZERO-I1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ), WORK( IW ), 1 );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        dswap(I2-IZERO+1, A( IOFF ), 1, WORK( IW ), 1 );
                     }
                  }

                  // Do for each value of NB in NBVAL

                  for (INB = 1; INB <= NNB; INB++) { // 50
                     NB = NBVAL( INB );
                     xlaenv(1, NB );

                     // Compute the L*L' or U'*U factorization of the band
                     // matrix.

                     dlacpy('Full', KD+1, N, A, LDAB, AFAC, LDAB );
                     SRNAMT = 'DPBTRF';
                     dpbtrf(UPLO, N, KD, AFAC, LDAB, INFO );

                     // Check error code from DPBTRF.

                     if ( INFO != IZERO ) {
                        alaerh(PATH, 'DPBTRF', INFO, IZERO, UPLO, N, N, KD, KD, NB, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 50;
                     }

                     // Skip the tests if INFO is not 0.

                     if (INFO != 0) GO TO 50;

*+    TEST 1
                     // Reconstruct matrix from factors and compute
                     // residual.

                     dlacpy('Full', KD+1, N, AFAC, LDAB, AINV, LDAB );
                     dpbt01(UPLO, N, KD, A, LDAB, AINV, LDAB, RWORK, RESULT( 1 ) );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 1 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, KD, NB, IMAT, 1, RESULT( 1 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;

                     // Only do other tests if this is the first blocksize.

                     if (INB > 1) GO TO 50;

                     // Form the inverse of A so we can get a good estimate
                     // of RCONDC = 1/(norm(A) * norm(inv(A))).

                     dlaset('Full', N, N, ZERO, ONE, AINV, LDA );
                     SRNAMT = 'DPBTRS';
                     dpbtrs(UPLO, N, KD, N, AFAC, LDAB, AINV, LDA, INFO );

                     // Compute RCONDC = 1/(norm(A) * norm(inv(A))).

                     ANORM = DLANSB( '1', UPLO, N, KD, A, LDAB, RWORK );
                     AINVNM = DLANGE( '1', N, N, AINV, LDA, RWORK );
                     if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                        RCONDC = ONE;
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM;
                     }

                     for (IRHS = 1; IRHS <= NNS; IRHS++) { // 40
                        NRHS = NSVAL( IRHS );

*+    TEST 2
                     // Solve and compute residual for A * X = B.

                        SRNAMT = 'DLARHS';
                        dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A, LDAB, XACT, LDA, B, LDA, ISEED, INFO );
                        dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'DPBTRS';
                        dpbtrs(UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA, INFO );

                     // Check error code from DPBTRS.

                        if (INFO != 0) CALL ALAERH( PATH, 'DPBTRS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        dpbt02(UPLO, N, KD, NRHS, A, LDAB, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

*+    TEST 3
                     // Check solution from generated exact solution.

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, and 6
                     // Use iterative refinement to improve the solution.

                        SRNAMT = 'DPBRFS';
                        dpbrfs(UPLO, N, KD, NRHS, A, LDAB, AFAC, LDAB, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO );

                     // Check error code from DPBRFS.

                        if (INFO != 0) CALL ALAERH( PATH, 'DPBRFS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                        dpbt05(UPLO, N, KD, NRHS, A, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 2; K <= 6; K++) { // 30
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, N, KD, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 30
                        NRUN = NRUN + 5;
                     } // 40

*+    TEST 7
                     // Get an estimate of RCOND = 1/CNDNUM.

                     SRNAMT = 'DPBCON';
                     dpbcon(UPLO, N, KD, AFAC, LDAB, ANORM, RCOND, WORK, IWORK, INFO );

                     // Check error code from DPBCON.

                     if (INFO != 0) CALL ALAERH( PATH, 'DPBCON', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );

                     RESULT( 7 ) = DGET06( RCOND, RCONDC );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 7 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )UPLO, N, KD, IMAT, 7, RESULT( 7 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;
                  } // 50
               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NB=', I4, ', type ', I2, ', test ', I2, ', ratio= ', G12.5 );
 9998 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') = ', G12.5 );
 9997 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ',', 10X, ' type ', I2, ', test(', I2, ') = ', G12.5 );
      RETURN;

      // End of DCHKPB

      }
