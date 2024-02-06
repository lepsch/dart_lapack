
// > \par Purpose:
// =============
// >
// > \verbatim
// >
// > SCHKQRT tests SGEQRT and SGEMQRT.
// > \endverbatim

// Arguments:
// ==========

// > \param[in] THRESH
// > \verbatim
// >          THRESH is REAL
// >          The threshold value for the test ratios.  A result is
// >          included in the output file if RESULT >= THRESH.  To have
// >          every test ratio printed, use THRESH = 0.
// > \endverbatim
// >
// > \param[in] TSTERR
// > \verbatim
// >          TSTERR is LOGICAL
// >          Flag that indicates whether error exits are to be tested.
// > \endverbatim
// >
// > \param[in] NM
// > \verbatim
// >          NM is INTEGER
// >          The number of values of M contained in the vector MVAL.
// > \endverbatim
// >
// > \param[in] MVAL
// > \verbatim
// >          MVAL is INTEGER array, dimension (NM)
// >          The values of the matrix row dimension M.
// > \endverbatim
// >
// > \param[in] NN
// > \verbatim
// >          NN is INTEGER
// >          The number of values of N contained in the vector NVAL.
// > \endverbatim
// >
// > \param[in] NVAL
// > \verbatim
// >          NVAL is INTEGER array, dimension (NN)
// >          The values of the matrix column dimension N.
// > \endverbatim
// >
// > \param[in] NNB
// > \verbatim
// >          NNB is INTEGER
// >          The number of values of NB contained in the vector NBVAL.
// > \endverbatim
// >
// > \param[in] NBVAL
// > \verbatim
// >          NBVAL is INTEGER array, dimension (NNB)
// >          The values of the blocksize NB.
// > \endverbatim
// >
// > \param[in] NOUT
// > \verbatim
// >          NOUT is INTEGER
// >          The unit number for output.
// > \endverbatim

// Authors:
// ========

// > \author Univ. of Tennessee
// > \author Univ. of California Berkeley
// > \author Univ. of Colorado Denver
// > \author NAG Ltd.

// > \ingroup single_lin

// =====================================================================
      void schkqrt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double               THRESH;
      int                MVAL( * ), NBVAL( * ), NVAL( * );

      int                NTESTS;
      const              NTESTS = 6 ;
      String             PATH;
      int                I, J, K, T, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRQRT, SQRT04
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Initialize constants

      PATH[1: 1] = 'S';
      PATH[2: 3] = 'QT';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      if (TSTERR) serrqrt( PATH, NOUT );
      INFOT = 0;

      // Do for each value of M in MVAL.

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );

         // Do for each possible value of NB

            MINMN = min( M, N );
            for (K = 1; K <= NNB; K++) {
               NB = NBVAL( K );
               if ( (NB <= MINMN) && (NB > 0) ) {

                  // Test SGEQRT and SGEMQRT

                  sqrt04(M, N, NB, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                  for (T = 1; T <= NTESTS; T++) {
                     if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )M, N, NB, T, RESULT( T );
                        NFAIL = NFAIL + 1;
                     }
                  }
                  NRUN = NRUN + NTESTS;
               }
            }
         }
      }

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=${.i5}, N=${.i5}, NB=${.i4} test(${.i2})=${.g12_5};
      return;
      }
