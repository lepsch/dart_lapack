      void zchkq3(final Array<bool> DOTYPE_, final int NM, final Array<int> MVAL_, final int NN, final Array<int> NVAL_, final int NNB, final Array<int> NBVAL_, final Array<int> NXVAL_, final int THRESH, final Array<double> A_, final Array<double> COPYA_, final Array<double> S_, final Array<double> TAU_, final Array<double> WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      double             S( * ), RWORK( * );
      Complex         A( * ), COPYA( * ), TAU( * ), WORK( * );
      // ..

      int                NTYPES;
      const              NTYPES = 6 ;
      int                NTESTS;
      const              NTESTS = 3 ;
      double             ONE, ZERO;
      Complex         CZERO;
      const              ONE = 1.0, ZERO = 0.0, CZERO = ( 0.0, 0.0 ) ;
      String             PATH;
      int                I, IHIGH, ILOW, IM, IMODE, IN, INB, INFO, ISTEP, K, LDA, LW, LWORK, M, MNMIN, MODE, N, NB, NRUN, NX;
      double             EPS;
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZQPT01, ZQRT11, ZQRT12;
      // EXTERNAL DLAMCH, ZQPT01, ZQRT11, ZQRT12
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, DLAORD, ICOPY, XLAENV, ZGEQP3, ZLACPY, ZLASET, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, IOUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];

      // Initialize constants and the random number seed.

      final PATH = '${'Zomplex precision'[0]}Q3';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10
      EPS = dlamch( 'Epsilon' );
      infoc.INFOT = 0;

      for (IM = 1; IM <= NM; IM++) { // 90

         // Do for each value of M in MVAL.

         final M = MVAL[IM];
         final LDA = max( 1, M );

         for (IN = 1; IN <= NN; IN++) { // 80

            // Do for each value of N in NVAL.

            final N = NVAL[IN];
            MNMIN = min( M, N );
            LWORK = max( 1, M*max( M, N )+4*MNMIN+max( M, N ) );

            for (IMODE = 1; IMODE <= NTYPES; IMODE++) { // 70
               if( !DOTYPE( IMODE ) ) GO TO 70;

               // Do for each type of matrix
               //    1:  zero matrix
               //    2:  one small singular value
               //    3:  geometric distribution of singular values
               //    4:  first n/2 columns fixed
               //    5:  last n/2 columns fixed
               //    6:  every second column fixed

               MODE = IMODE;
               if (IMODE > 3) MODE = 1;

               // Generate test matrix of size m by n using
               // singular value distribution indicated by `mode'.

               for (I = 1; I <= N; I++) { // 20
                  IWORK[I] = 0;
               } // 20
               if ( IMODE == 1 ) {
                  zlaset('Full', M, N, CZERO, CZERO, COPYA, LDA );
                  for (I = 1; I <= MNMIN; I++) { // 30
                     S[I] = ZERO;
                  } // 30
               } else {
                  zlatms(M, N, 'Uniform', ISEED, 'Nonsymm', S, MODE, ONE / EPS, ONE, M, N, 'No packing', COPYA, LDA, WORK, INFO );
                  if ( IMODE >= 4 ) {
                     if ( IMODE == 4 ) {
                        ILOW = 1;
                        ISTEP = 1;
                        IHIGH = max( 1, N / 2 );
                     } else if ( IMODE == 5 ) {
                        ILOW = max( 1, N / 2 );
                        ISTEP = 1;
                        IHIGH = N;
                     } else if ( IMODE == 6 ) {
                        ILOW = 1;
                        ISTEP = 2;
                        IHIGH = N;
                     }
                     for (I = ILOW; ISTEP < 0 ? I >= IHIGH : I <= IHIGH; I += ISTEP) { // 40
                        IWORK[I] = 1;
                     } // 40
                  }
                  dlaord('Decreasing', MNMIN, S, 1 );
               }

               for (INB = 1; INB <= NNB; INB++) { // 60

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  final NB = NBVAL[INB];
                  xlaenv(1, NB );
                  NX = NXVAL( INB );
                  xlaenv(3, NX );

                  // Save A and its singular values and a copy of
                  // vector IWORK.

                  zlacpy('All', M, N, COPYA, LDA, A, LDA );
                  icopy(N, IWORK( 1 ), 1, IWORK( N+1 ), 1 );

                  // Workspace needed.

                  LW = NB*( N+1 );

                 srnamc.SRNAMT = 'ZGEQP3';
                  zgeqp3(M, N, A, LDA, IWORK( N+1 ), TAU, WORK, LW, RWORK, INFO );

                  // Compute norm(svd(a) - svd(r))

                  RESULT[1] = ZQRT12( M, N, A, LDA, S, WORK, LWORK, RWORK );

                  // Compute norm( A*P - Q*R )

                  RESULT[2] = ZQPT01( M, N, MNMIN, COPYA, A, LDA, TAU, IWORK( N+1 ), WORK, LWORK );

                  // Compute Q'*Q

                  RESULT[3] = ZQRT11( M, MNMIN, A, LDA, TAU, WORK, LWORK );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NTESTS; K++) { // 50
                     if ( RESULT[K] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                        NOUT.println( 9999 )'ZGEQP3', M, N, NB, IMODE, K, RESULT( K );
                        NFAIL++;
                     }
                  } // 50
                  NRUN +=  NTESTS;

               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${} M =${M.i5}, N =${N.i5}, NB =${NB.i4}, type ${IMAT.i2}, test ${.i2}, ratio =${RESULT[].g12_5};
      }
