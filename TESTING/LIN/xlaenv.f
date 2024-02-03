      SUBROUTINE XLAENV( ISPEC, NVALUE )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                ISPEC, NVALUE;
      // ..
*
*  =====================================================================
*
      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      COMMON             / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /
      // ..
      // .. Executable Statements ..
*
      IF( ISPEC.GE.1 .AND. ISPEC.LE.9 ) THEN
         IPARMS( ISPEC ) = NVALUE
      END IF
*
      RETURN
*
      // End of XLAENV
*
      END
