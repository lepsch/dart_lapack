      SUBROUTINE ZSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL TRYRAC
*     ..
*     .. External Subroutines ..
      EXTERNAL ZSTEMR
*     ..
*     .. Executable Statements ..
      INFO = 0
      TRYRAC = .FALSE.
       CALL ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )
*
*     End of ZSTEGR
*
      END
