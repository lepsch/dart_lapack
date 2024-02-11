import 'common.dart';
      void xlaenv(final int ISPEC, final int NVALUE,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

      if ( ISPEC >= 1 && ISPEC <= 9 ) {
         claenv.IPARMS[ISPEC] = NVALUE;
      }

      }
