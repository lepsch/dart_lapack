      SUBROUTINE CSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      REAL             ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      REAL               D( * ), E( * ), W( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      bool    TRYRAC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSTEMR
      // ..
      // .. Executable Statements ..
      INFO = 0
      TRYRAC = .FALSE.
       CALL CSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )
*
      // End of CSTEGR
*
      END
