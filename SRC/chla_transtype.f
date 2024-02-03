      String      FUNCTION CHLA_TRANSTYPE( TRANS );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                TRANS;
      // ..

*  =====================================================================

      // .. Parameters ..
      int     BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS;
      const     BLAS_NO_TRANS = 111, BLAS_TRANS = 112, BLAS_CONJ_TRANS = 113 ;
      // ..
      // .. Executable Statements ..
      IF( TRANS.EQ.BLAS_NO_TRANS ) THEN
         CHLA_TRANSTYPE = 'N'
      ELSE IF( TRANS.EQ.BLAS_TRANS ) THEN
         CHLA_TRANSTYPE = 'T'
      ELSE IF( TRANS.EQ.BLAS_CONJ_TRANS ) THEN
         CHLA_TRANSTYPE = 'C'
      ELSE
         CHLA_TRANSTYPE = 'X'
      END IF
      RETURN

      // End of CHLA_TRANSTYPE

      }
