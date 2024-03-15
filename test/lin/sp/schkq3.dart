      void schkq3(final int DOTYPE, final int NM, final int MVAL, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NXVAL, final int THRESH, final int A, final int COPYA, final int S, final int TAU, final Array<double> _WORK_, final Array<int> IWORK_, final int NOUT,) {
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NM, NN, NNB, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      double               A( * ), COPYA( * ), S( * ), TAU( * ), WORK( * );
      // ..

      int                NTYPES;
      const              NTYPES = 6 ;
      int                NTESTS;
      const              NTESTS = 3 ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             PATH;
      int                I, IHIGH, ILOW, IM, IMODE, IN, INB, INFO, ISTEP, K, LDA, LW, LWORK, M, MNMIN, MODE, N, NB, NERRS, NFAIL, NRUN, NX;
      double               EPS;
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SQPT01, SQRT11, SQRT12;
      // EXTERNAL SLAMCH, SQPT01, SQRT11, SQRT12
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, ICOPY, SGEQP3, SLACPY, SLAORD, SLASET, SLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, IOUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'Q3';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10
      EPS = SLAMCH( 'Epsilon' );
      INFOT = 0;

      for (IM = 1; IM <= NM; IM++) { // 90

         // Do for each value of M in MVAL.

         M = MVAL( IM );
         LDA = max( 1, M );

         for (IN = 1; IN <= NN; IN++) { // 80

            // Do for each value of N in NVAL.

            N = NVAL( IN );
            MNMIN = min( M, N );
            LWORK = max( 1, M*max( M, N )+4*MNMIN+max( M, N ), M*N + 2*MNMIN + 4*N );

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
                  slaset('Full', M, N, ZERO, ZERO, COPYA, LDA );
                  for (I = 1; I <= MNMIN; I++) { // 30
                     S[I] = ZERO;
                  } // 30
               } else {
                  slatms(M, N, 'Uniform', ISEED, 'Nonsymm', S, MODE, ONE / EPS, ONE, M, N, 'No packing', COPYA, LDA, WORK, INFO );
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
                  slaord('Decreasing', MNMIN, S, 1 );
               }

               for (INB = 1; INB <= NNB; INB++) { // 60

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );
                  NX = NXVAL( INB );
                  xlaenv(3, NX );

                  // Get a working copy of COPYA into A and a copy of
                  // vector IWORK.

                  slacpy('All', M, N, COPYA, LDA, A, LDA );
                  icopy(N, IWORK( 1 ), 1, IWORK( N+1 ), 1 );

                  // Compute the QR factorization with pivoting of A

                  LW = max( 1, 2*N+NB*( N+1 ) );

                  // Compute the QP3 factorization of A

                 srnamc.SRNAMT = 'SGEQP3';
                  sgeqp3(M, N, A, LDA, IWORK( N+1 ), TAU, WORK, LW, INFO );

                  // Compute norm(svd(a) - svd(r))

                  RESULT[1] = SQRT12( M, N, A, LDA, S, WORK, LWORK );

                  // Compute norm( A*P - Q*R )

                  RESULT[2] = SQPT01( M, N, MNMIN, COPYA, A, LDA, TAU, IWORK( N+1 ), WORK, LWORK );

                  // Compute Q'*Q

                  RESULT[3] = SQRT11( M, MNMIN, A, LDA, TAU, WORK, LWORK );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NTESTS; K++) { // 50
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )'SGEQP3', M, N, NB, IMODE, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 50
                  NRUN = NRUN + NTESTS;

               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${} M =${.i5}, N =${.i5}, NB =${.i4}, type ${.i2}, test ${.i2}, ratio =${.g12_5};
      }