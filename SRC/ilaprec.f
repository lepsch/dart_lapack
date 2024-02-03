      int     FUNCTION ILAPREC( PREC );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PREC;
      // ..

*  =====================================================================

      // .. Parameters ..
      int     BLAS_PREC_SINGLE, BLAS_PREC_DOUBLE, BLAS_PREC_INDIGENOUS, BLAS_PREC_EXTRA;
      const     BLAS_PREC_SINGLE = 211, BLAS_PREC_DOUBLE = 212, BLAS_PREC_INDIGENOUS = 213, BLAS_PREC_EXTRA = 214 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
      if ( LSAME( PREC, 'S' ) ) {
         ILAPREC = BLAS_PREC_SINGLE;
      } else if ( LSAME( PREC, 'D' ) ) {
         ILAPREC = BLAS_PREC_DOUBLE;
      } else if ( LSAME( PREC, 'I' ) ) {
         ILAPREC = BLAS_PREC_INDIGENOUS;
      } else if ( LSAME( PREC, 'X' ) || LSAME( PREC, 'E' ) ) {
         ILAPREC = BLAS_PREC_EXTRA;
      } else {
         ILAPREC = -1;
      }
      return;

      // End of ILAPREC

      }
