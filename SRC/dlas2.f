      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             F, G, H, SSMAX, SSMIN;
      // ..

*  ====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             ONE;
      const              ONE = 1.0D0 ;
      double             TWO;
      const              TWO = 2.0D0 ;
      // ..
      // .. Local Scalars ..
      double             AS, AT, AU, C, FA, FHMN, FHMX, GA, HA;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      FA = ABS( F )
      GA = ABS( G )
      HA = ABS( H )
      FHMN = MIN( FA, HA )
      FHMX = MAX( FA, HA )
      if ( FHMN == ZERO ) {
         SSMIN = ZERO
         if ( FHMX == ZERO ) {
            SSMAX = GA
         } else {
            SSMAX = MAX( FHMX, GA )*SQRT( ONE+ ( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
         }
      } else {
         if ( GA.LT.FHMX ) {
            AS = ONE + FHMN / FHMX
            AT = ( FHMX-FHMN ) / FHMX
            AU = ( GA / FHMX )**2
            C = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
            SSMIN = FHMN*C
            SSMAX = FHMX / C
         } else {
            AU = FHMX / GA
            if ( AU == ZERO ) {

               // Avoid possible harmful underflow if exponent range
               // asymmetric (true SSMIN may not underflow even if
               // AU underflows)

               SSMIN = ( FHMN*FHMX ) / GA
               SSMAX = GA
            } else {
               AS = ONE + FHMN / FHMX
               AT = ( FHMX-FHMN ) / FHMX
               C = ONE / ( SQRT( ONE+( AS*AU )**2 )+ SQRT( ONE+( AT*AU )**2 ) )
               SSMIN = ( FHMN*C )*AU
               SSMIN = SSMIN + SSMIN
               SSMAX = GA / ( C+C )
            }
         }
      }
      RETURN

      // End of DLAS2

      }
