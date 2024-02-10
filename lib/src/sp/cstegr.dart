      void cstegr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, final Matrix<double> Z, final int LDZ, ISUPPZ, final Array<double> WORK, final int LWORK, final Array<int> IWORK, final int LIWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      int                ISUPPZ( * ), IWORK( * );
      double               D( * ), E( * ), W( * ), WORK( * );
      Complex            Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool    TRYRAC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSTEMR
      INFO = 0;
      TRYRAC = false;
       cstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO );
      }
