      void zchklqtp(final double THRESH, final bool TSTERR, final int NM, final Array<int> MVAL_, final int NN, final Array<int> NVAL_, final int NNB, final Array<int> NBVAL_, final Nout NOUT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

      int                NTESTS;
      const              NTESTS = 6 ;
      String             PATH;
      int                I, J, K, L, T, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZERRLQTP, ZLQT04
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Initialize constants

      PATH[1: 1] = 'Z';
      PATH[2: 3] = 'XQ';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);

      // Test the error exits

      if (TSTERR) zerrlqtp( PATH, NOUT );
      infoc.INFOT = 0;

      // Do for each value of M

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );

            // Do for each value of L

            MINMN = min( M, N );
            for (L = 0; max( MINMN, 1 ) < 0 ? L >= MINMN : L <= MINMN; L += max( MINMN, 1 )) {

               // Do for each possible value of NB

               for (K = 1; K <= NNB; K++) {
                  NB = NBVAL( K );

                  // Test DTPLQT and DTPMLQT

                  if ( (NB <= M) && (NB > 0) ) {
                     zlqt05(M, N, L, NB, RESULT );

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (T = 1; T <= NTESTS; T++) {
                     if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                           NOUT.println( 9999 )M, N, NB, L, T, RESULT( T );
                           NFAIL++;
                        }
                     }
                     NRUN +=  NTESTS;
                  }
               }
            }
         }
      }

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=${M.i5}, N=${N.i5}, NB=', I4,' L=${.i4} test(${.i2})=${.g12_5};
      }
