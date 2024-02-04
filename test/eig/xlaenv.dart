      // int                IPARMS( 100 );

      void xlaenv(final int ISPEC, final int NVALUE, ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /;

      if ( ISPEC >= 1 && ISPEC <= 16 ) {
         IPARMS[ISPEC] = NVALUE;
      }

      return;
      }
