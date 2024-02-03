      SUBROUTINE XLAENV( ISPEC, NVALUE )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ISPEC, NVALUE;
      // ..

*  =====================================================================

      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /
      // ..
      // .. Executable Statements ..

      if ( ISPEC >= 1 && ISPEC.LE.9 ) {
         IPARMS( ISPEC ) = NVALUE
      }

      RETURN

      // End of XLAENV

      }
