      SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, LDC, M, N;
      COMPLEX            TAU
      // ..
      // .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      const              ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               APPLYLEFT;
      int                I, LASTV, LASTC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERC
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILACLR, ILACLC;
      // EXTERNAL LSAME, ILACLR, ILACLC
      // ..
      // .. Executable Statements ..

      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      if ( TAU.NE.ZERO ) {
      // Set up variables for scanning V.  LASTV begins pointing to the end
      // of V.
         if ( APPLYLEFT ) {
            LASTV = M
         } else {
            LASTV = N
         }
         if ( INCV.GT.0 ) {
            I = 1 + (LASTV-1) * INCV
         } else {
            I = 1
         }
      // Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         if ( APPLYLEFT ) {
      // Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILACLC(LASTV, N, C, LDC)
         } else {
      // Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILACLR(M, LASTV, C, LDC)
         }
      }
      // Note that lastc.eq.0 renders the BLAS operations null; no special
      // case is needed at this level.
      if ( APPLYLEFT ) {

         // Form  H * C

         if ( LASTV.GT.0 ) {

            // w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)

            cgemv('Conjugate transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

            // C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H

            cgerc(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC );
         }
      } else {

         // Form  C * H

         if ( LASTV.GT.0 ) {

            // w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)

            cgemv('No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

            // C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H

            cgerc(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC );
         }
      }
      RETURN

      // End of CLARF

      }
