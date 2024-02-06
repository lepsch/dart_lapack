      void zchkps(DOTYPE, NN, NVAL, NNB, NBVAL, NRANK, RANKVAL, THRESH, TSTERR, NMAX, A, AFAC, PERM, PIV, WORK, RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double             THRESH;
      int                NMAX, NN, NNB, NOUT, NRANK;
      bool               TSTERR;
      Complex         A( * ), AFAC( * ), PERM( * ), WORK( * );
      double             RWORK( * );
      int                NBVAL( * ), NVAL( * ), PIV( * ), RANKVAL( * );
      bool               DOTYPE( * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      double             ANORM, CNDNUM, RESULT, TOL;
      int                COMPRANK, I, IMAT, IN, INB, INFO, IRANK, IUPLO, IZERO, KL, KU, LDA, MODE, N, NB, NERRS, NFAIL, NIMAT, NRUN, RANK, RANKDIFF;
      String             DIST, TYPE, UPLO;
      String             PATH;
      int                ISEED( 4 ), ISEEDY( 4 );
      String             UPLOS( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRPS, ZLACPY, ZLATB5, ZLATMT, ZPST01, ZPSTRF
      // ..
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NUNIT;
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, CEILING
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex Precision';
      PATH[2: 3] = 'PS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 100
         ISEED[I] = ISEEDY( I );
      } // 100

      // Test the error exits

      if (TSTERR) zerrps( PATH, NOUT );
      infoc.INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 150
         N = NVAL( IN );
         LDA = max( N, 1 );
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;
         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 140

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 140;

               // Do for each value of RANK in RANKVAL

            for (IRANK = 1; IRANK <= NRANK; IRANK++) { // 130

               // Only repeat test 3 to 5 for different ranks
               // Other tests use full rank

               if( ( IMAT < 3 || IMAT > 5 ) && IRANK > 1 ) GO TO 130;

               RANK = CEILING( ( N * (RANKVAL( IRANK )).toDouble() ) / 100.0 );


            // Do first for UPLO = 'U', then for UPLO = 'L'

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 120
                  UPLO = UPLOS( IUPLO );

               // Set up parameters with ZLATB5 and generate a test matrix
               // with ZLATMT.

                  zlatb5(PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                 srnamc.SRNAMT = 'ZLATMT';
                  zlatmt(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, RANK, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from ZLATMT.

                  if ( INFO != 0 ) {
                    alaerh(PATH, 'ZLATMT', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 120;
                  }

               // Do for each value of NB in NBVAL

                  for (INB = 1; INB <= NNB; INB++) { // 110
                     NB = NBVAL( INB );
                     xlaenv(1, NB );

                  // Compute the pivoted L*L' or U'*U factorization
                  // of the matrix.

                     zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                    srnamc.SRNAMT = 'ZPSTRF';

                  // Use default tolerance

                     TOL = -ONE;
                     zpstrf(UPLO, N, AFAC, LDA, PIV, COMPRANK, TOL, RWORK, INFO );

                  // Check error code from ZPSTRF.

                     if( (INFO < IZERO) || (INFO != IZERO && RANK == N) || (INFO <= IZERO && RANK < N) ) THEN;
                        alaerh(PATH, 'ZPSTRF', INFO, IZERO, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 110;
                     }

                  // Skip the test if INFO is not 0.

                     if (INFO != 0) GO TO 110;

                  // Reconstruct matrix from factors and compute residual.

                  // PERM holds permuted L*L^T or U^T*U

                     zpst01(UPLO, N, A, LDA, AFAC, LDA, PERM, LDA, PIV, RWORK, RESULT, COMPRANK );

                  // Print information about the tests that did not pass
                  // the threshold or where computed rank was not RANK.

                     if (N == 0) COMPRANK = 0;
                     RANKDIFF = RANK - COMPRANK;
                     if ( RESULT >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )UPLO, N, RANK, RANKDIFF, NB, IMAT, RESULT;
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;
                  } // 110

               } // 120
            } // 130
         } // 140
      } // 150

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', RANK =', I3, ', Diff =', I5, ', NB =', I4, ', type ', I2, ', Ratio =', G12.5 );
      return;
      }
