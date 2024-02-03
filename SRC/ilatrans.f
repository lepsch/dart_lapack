      int     FUNCTION ILATRANS( TRANS );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      // ..

*  =====================================================================

      // .. Parameters ..
      int     BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS;
      PARAMETER ( BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
      IF( LSAME( TRANS, 'N' ) ) THEN
         ILATRANS = BLAS_NO_TRANS
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         ILATRANS = BLAS_TRANS
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
         ILATRANS = BLAS_CONJ_TRANS
      ELSE
         ILATRANS = -1
      END IF
      RETURN

      // End of ILATRANS

      }
