      void cchklq(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC, B, X, XACT, TAU, WORK, RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NOUT, NRHS;
      double               THRESH;
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      double               RWORK( * );
      Complex            A( * ), AC( * ), AF( * ), AL( * ), AQ( * ), B( * ), TAU( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      String             DIST, TYPE;
      String             PATH;
      int                I, IK, IM, IMAT, IN, INB, INFO, K, KL, KU, LDA, LWORK, M, MINMN, MODE, N, NB, NERRS, NFAIL, NK, NRUN, NT, NX;
      double               ANORM, CNDNUM;
      int                ISEED( 4 ), ISEEDY( 4 ), KVAL( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CERRLQ, CGELS, CGET02, CLACPY, CLARHS, CLATB4, CLATMS, CLQT01, CLQT02, CLQT03, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'LQ';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrlq( PATH, NOUT );
      INFOT = 0;
      xlaenv(2, 2 );

      LDA = NMAX;
      LWORK = NMAX*max( NMAX, NRHS );

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 70
         M = MVAL( IM );

         // Do for each value of N in NVAL.

         for (IN = 1; IN <= NN; IN++) { // 60
            N = NVAL( IN );
            MINMN = min( M, N );
            for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 50

               // Do the tests only if DOTYPE( IMAT ) is true.

               if( !DOTYPE( IMAT ) ) GO TO 50;

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

              srnamc.SRNAMT = 'CLATMS';
               clatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 50;
               }

               // Set some values for K: the first value must be MINMN,
               // corresponding to the call of CLQT01; other values are
               // used in the calls of CLQT02, and must not exceed MINMN.

               KVAL[1] = MINMN;
               KVAL[2] = 0;
               KVAL[3] = 1;
               KVAL[4] = MINMN / 2;
               if ( MINMN == 0 ) {
                  NK = 1;
               } else if ( MINMN == 1 ) {
                  NK = 2;
               } else if ( MINMN <= 3 ) {
                  NK = 3;
               } else {
                  NK = 4;
               }

               // Do for each value of K in KVAL

               for (IK = 1; IK <= NK; IK++) { // 40
                  K = KVAL( IK );

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  for (INB = 1; INB <= NNB; INB++) { // 30
                     NB = NBVAL( INB );
                     xlaenv(1, NB );
                     NX = NXVAL( INB );
                     xlaenv(3, NX );
                     for (I = 1; I <= NTESTS; I++) {
                        RESULT[I] = ZERO;
                     }
                     NT = 2;
                     if ( IK == 1 ) {

                        // Test CGELQF

                        clqt01(M, N, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );
                     } else if ( M <= N ) {

                        // Test CUNGLQ, using factorization
                        // returned by CLQT01

                        clqt02(M, N, K, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );
                     }
                     if ( M >= K ) {

                        // Test CUNMLQ, using factorization returned
                        // by CLQT01

                        clqt03(M, N, K, AF, AC, AL, AQ, LDA, TAU, WORK, LWORK, RWORK, RESULT( 3 ) );
                        NT = NT + 4;

                        // If M<=N and K=M, call CGELS to solve a system
                        // with NRHS right hand sides and compute the
                        // residual.

                        if ( K == M && INB == 1 ) {

                           // Generate a solution and set the right
                           // hand side.

                          srnamc.SRNAMT = 'CLARHS';
                           clarhs(PATH, 'New', 'Full', 'No transpose', M, N, 0, 0, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );

                           clacpy('Full', M, NRHS, B, LDA, X, LDA );

                           // Reset AF to the original matrix. CGELS
                           // factors the matrix before solving the system.

                           clacpy('Full', M, N, A, LDA, AF, LDA );

                          srnamc.SRNAMT = 'CGELS';
                           cgels('No transpose', M, N, NRHS, AF, LDA, X, LDA, WORK, LWORK, INFO );

                           // Check error code from CGELS.

                           if (INFO != 0) alaerh( PATH, 'CGELS', INFO, 0, 'N', M, N, NRHS, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                           cget02('No transpose', M, N, NRHS, A, LDA, X, LDA, B, LDA, RWORK, RESULT( 7 ) );
                           NT = NT + 1;
                        }
                     }

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (I = 1; I <= NT; I++) { // 20
                        if ( RESULT( I ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9999 )M, N, K, NB, NX, IMAT, I, RESULT( I );
                           NFAIL = NFAIL + 1;
                        }
                     } // 20
                     NRUN = NRUN + NT;
                  } // 30
               } // 40
            } // 50
         } // 60
      } // 70

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=${.i5}, N=${.i5}, K=${.i5}, NB=${.i4}, NX=${.i5}, type ${.i2}, test(${.i2})=${.g12_5};
      return;
      }
