import 'common.dart';
      void xlaenv(ISPEC, NVALUE ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ISPEC, NVALUE;
      // ..

// =====================================================================

      // .. Arrays in Common ..
      int                claenv.IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / claenv.IPARMS
      // ..
      // .. Save statement ..
      // SAVE               / CLAENV /;
      // ..
      // .. Executable Statements ..

      if ( ISPEC >= 1 && ISPEC <= 9 ) {
         claenv.IPARMS[ISPEC] = NVALUE;
      }

      return;
      }
