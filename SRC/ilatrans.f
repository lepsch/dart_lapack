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
      const     BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
      if ( LSAME( TRANS, 'N' ) ) {
         ILATRANS = BLAS_NO_TRANS;
      } else if ( LSAME( TRANS, 'T' ) ) {
         ILATRANS = BLAS_TRANS;
      } else if ( LSAME( TRANS, 'C' ) ) {
         ILATRANS = BLAS_CONJ_TRANS;
      } else {
         ILATRANS = -1;
      }
      RETURN;

      // End of ILATRANS

      }
