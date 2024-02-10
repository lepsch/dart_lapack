      void slas2(F, G, H, SSMIN, final int SSMAX) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               F, G, H, SSMAX, SSMIN;
      // ..

// ====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      double               TWO;
      const              TWO = 2.0 ;
      double               AS, AT, AU, C, FA, FHMN, FHMX, GA, HA;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT

      FA = ( F ).abs();
      GA = ( G ).abs();
      HA = ( H ).abs();
      FHMN = min( FA, HA );
      FHMX = max( FA, HA );
      if ( FHMN == ZERO ) {
         SSMIN = ZERO;
         if ( FHMX == ZERO ) {
            SSMAX = GA;
         } else {
            SSMAX = max( FHMX, GA )*sqrt( ONE+ ( min( FHMX, GA ) / max( FHMX, GA ) )**2 );
         }
      } else {
         if ( GA < FHMX ) {
            AS = ONE + FHMN / FHMX;
            AT = ( FHMX-FHMN ) / FHMX;
            AU = ( GA / FHMX )**2;
            C = TWO / ( sqrt( AS*AS+AU )+sqrt( AT*AT+AU ) );
            SSMIN = FHMN*C;
            SSMAX = FHMX / C;
         } else {
            AU = FHMX / GA;
            if ( AU == ZERO ) {

               // Avoid possible harmful underflow if exponent range
               // asymmetric (true SSMIN may not underflow even if
               // AU underflows)

               SSMIN = ( FHMN*FHMX ) / GA;
               SSMAX = GA;
            } else {
               AS = ONE + FHMN / FHMX;
               AT = ( FHMX-FHMN ) / FHMX;
               C = ONE / ( sqrt( ONE+( AS*AU )**2 )+ sqrt( ONE+( AT*AU )**2 ) );
               SSMIN = ( FHMN*C )*AU;
               SSMIN = SSMIN + SSMIN;
               SSMAX = GA / ( C+C );
            }
         }
      }
      }
