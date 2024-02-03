      String      FUNCTION CHLA_TRANSTYPE( TRANS );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                TRANS;
      // ..

// =====================================================================

      // .. Parameters ..
      int     BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS;
      const     BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113 ;
      // ..
      // .. Executable Statements ..
      if ( TRANS == BLAS_NO_TRANS ) {
         CHLA_TRANSTYPE = 'N';
      } else if ( TRANS == BLAS_TRANS ) {
         CHLA_TRANSTYPE = 'T';
      } else if ( TRANS == BLAS_CONJ_TRANS ) {
         CHLA_TRANSTYPE = 'C';
      } else {
         CHLA_TRANSTYPE = 'X';
      }
      return;

      // End of CHLA_TRANSTYPE

      }
