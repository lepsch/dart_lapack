      int     FUNCTION ILAPREC( PREC );
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             PREC;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int     BLAS_PREC_SINGLE, BLAS_PREC_DOUBLE, BLAS_PREC_INDIGENOUS, BLAS_PREC_EXTRA       PARAMETER ( BLAS_PREC_SINGLE = 211, BLAS_PREC_DOUBLE = 212, BLAS_PREC_INDIGENOUS = 213, BLAS_PREC_EXTRA = 214 );
*     ..
*     .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( PREC, 'S' ) ) THEN
         ILAPREC = BLAS_PREC_SINGLE
      ELSE IF( LSAME( PREC, 'D' ) ) THEN
         ILAPREC = BLAS_PREC_DOUBLE
      ELSE IF( LSAME( PREC, 'I' ) ) THEN
         ILAPREC = BLAS_PREC_INDIGENOUS
      ELSE IF( LSAME( PREC, 'X' ) .OR. LSAME( PREC, 'E' ) ) THEN
         ILAPREC = BLAS_PREC_EXTRA
      ELSE
         ILAPREC = -1
      END IF
      RETURN
*
*     End of ILAPREC
*
      END
