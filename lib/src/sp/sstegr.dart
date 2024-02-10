      void sstegr(final int JOBZ, final int RANGE, final int N, final int D, final int E, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int M, final int W, final Matrix<double> Z, final int LDZ, final int ISUPPZ, final Array<double> WORK, final int LWORK, final Array<int> IWORK, final int LIWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      int                ISUPPZ( * ), IWORK( * );
      double               D( * ), E( * ), W( * ), WORK( * );
      double               Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool    TRYRAC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSTEMR
      INFO = 0;
      TRYRAC = false;
       sstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO );
      }
