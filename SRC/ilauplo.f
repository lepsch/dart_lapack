      int     FUNCTION ILAUPLO( UPLO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int     BLAS_UPPER, BLAS_LOWER
      PARAMETER ( BLAS_UPPER = 121, BLAS_LOWER = 122 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( UPLO, 'U' ) ) THEN
         ILAUPLO = BLAS_UPPER
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         ILAUPLO = BLAS_LOWER
      ELSE
         ILAUPLO = -1
      END IF
      RETURN
*
*     End of ILAUPLO
*
      END
