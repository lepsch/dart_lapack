      SUBROUTINE DSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      double           ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             D( * ), E( * ), W( * ), WORK( * );
      double             Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool    TRYRAC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSTEMR
      // ..
      // .. Executable Statements ..
      INFO = 0;
      TRYRAC = false;
       dstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO );

      // End of DSTEGR

      }
