import 'common.dart';

      void dchkqrtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

      const              NTESTS = 6 ;
      String             PATH;
      int                I, J, K, L, T, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRQRTP
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT

      // Initialize constants

      PATH[1: 1] = 'D';
      PATH[2: 3] = 'QX';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      if (TSTERR) derrqrtp( PATH, NOUT );
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

                  // Test DTPQRT and DTPMQRT

                  if ( (NB <= N) && (NB > 0) ) {
                     dqrt05(M, N, L, NB, RESULT );

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (T = 1; T <= NTESTS; T++) {
                     if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9999 )M, N, NB, L, T, RESULT( T );
                           NFAIL = NFAIL + 1;
                        }
                     }
                     NRUN = NRUN + NTESTS;
                  }
               }
            }
         }
      }

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=${.i5}, N=${.i5}, NB=', I4,' L=${.i4} test(${.i2})=${.g12_5};
      return;
      }
