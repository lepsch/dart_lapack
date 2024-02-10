      void schkps(DOTYPE, NN, NVAL, NNB, NBVAL, NRANK, RANKVAL, THRESH, TSTERR, NMAX, A, AFAC, PERM, PIV, final Array<double> _WORK, final Array<double> RWORK, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               THRESH;
      int                NMAX, NN, NNB, NOUT, NRANK;
      bool               TSTERR;
      double               A( * ), AFAC( * ), PERM( * ), RWORK( * ), WORK( * );
      int                NBVAL( * ), NVAL( * ), PIV( * ), RANKVAL( * );
      bool               DOTYPE( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      double               ANORM, CNDNUM, RESULT, TOL;
      int                COMPRANK, I, IMAT, IN, INB, INFO, IRANK, IUPLO, IZERO, KL, KU, LDA, MODE, N, NB, NERRS, NFAIL, NIMAT, NRUN, RANK, RANKDIFF;
      String             DIST, TYPE, UPLO;
      String             PATH;
      int                ISEED( 4 ), ISEEDY( 4 );
      String             UPLOS( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRPS, SLACPY, SLATB5, SLATMT, SPST01, SPSTRF, XLAENV
      // ..
      // .. Scalars in Common ..
      int                INFOT, NUNIT;
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, CEILING
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Single Precision';
      PATH[2: 3] = 'PS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 100
         ISEED[I] = ISEEDY( I );
      } // 100

      // Test the error exits

      if (TSTERR) serrps( PATH, NOUT );
      INFOT = 0;
      xlaenv(2, 2 );

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

               RANK = CEILING( ( N * double( RANKVAL( IRANK ) ) ) / 100.0 );


            // Do first for UPLO = 'U', then for UPLO = 'L'

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 120
                  UPLO = UPLOS( IUPLO );

               // Set up parameters with SLATB5 and generate a test matrix
               // with SLATMT.

                  slatb5(PATH, IMAT, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                 srnamc.SRNAMT = 'SLATMT';
                  slatmt(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, RANK, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from SLATMT.

                  if ( INFO != 0 ) {
                    alaerh(PATH, 'SLATMT', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 120;
                  }

               // Do for each value of NB in NBVAL

                  for (INB = 1; INB <= NNB; INB++) { // 110
                     NB = NBVAL( INB );
                     xlaenv(1, NB );

                  // Compute the pivoted L*L' or U'*U factorization
                  // of the matrix.

                     slacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                    srnamc.SRNAMT = 'SPSTRF';

                  // Use default tolerance

                     TOL = -ONE;
                     spstrf(UPLO, N, AFAC, LDA, PIV, COMPRANK, TOL, WORK, INFO );

                  // Check error code from SPSTRF.

                     if( (INFO < IZERO) || (INFO != IZERO && RANK == N) || (INFO <= IZERO && RANK < N) ) THEN;
                        alaerh(PATH, 'SPSTRF', INFO, IZERO, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 110;
                     }

                  // Skip the test if INFO is not 0.

                     if (INFO != 0) GO TO 110;

                  // Reconstruct matrix from factors and compute residual.

                  // PERM holds permuted L*L^T or U^T*U

                     spst01(UPLO, N, A, LDA, AFAC, LDA, PERM, LDA, PIV, RWORK, RESULT, COMPRANK );

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

 9999 FORMAT( ' UPLO = ''${.a1}'', N =${.i5}, RANK =${.i3}, Diff =${.i5}, NB =${.i4}, type ${.i2}, Ratio =${.g12_5};
      }
