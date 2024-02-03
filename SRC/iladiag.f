      int     FUNCTION ILADIAG( DIAG );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG;
      // ..

*  =====================================================================

      // .. Parameters ..
      int     BLAS_NON_UNIT_DIAG, BLAS_UNIT_DIAG;
      const     BLAS_NON_UNIT_DIAG = 131, BLAS_UNIT_DIAG = 132 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
      IF( LSAME( DIAG, 'N' ) ) THEN
         ILADIAG = BLAS_NON_UNIT_DIAG
      ELSE IF( LSAME( DIAG, 'U' ) ) THEN
         ILADIAG = BLAS_UNIT_DIAG
      } else {
         ILADIAG = -1
      END IF
      RETURN

      // End of ILADIAG

      }
