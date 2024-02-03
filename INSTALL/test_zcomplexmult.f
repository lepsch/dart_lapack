void main() {
*  -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // ..
      // .. Constants ..
      int               nNaN, nInf;
      const           nNaN = 3, nInf = 5 ;
      double complex    czero, cone
      const           czero = DCMPLX( 0.0d0, 0.0d0 ), cone  = DCMPLX( 1.0d0, 0.0d0 ) ;
      // ..
      // .. Local Variables ..
      int               i, nFailingTests, nTests;
      double            aInf, aNaN, OV;
      double complex    Y, R, cInf( nInf ), cNaN( nNaN )

      // .. Intrinsic Functions ..
      // intrinsic HUGE, DCMPLX


      // .. Initialize error counts ..
      nFailingTests = 0
      nTests = 0

      // .. Inf entries ..
      OV = HUGE(0.0d0)
      aInf = OV * 2
      cInf(1) = DCMPLX( aInf, 0.0d0 )
      cInf(2) = DCMPLX(-aInf, 0.0d0 )
      cInf(3) = DCMPLX( 0.0d0, aInf )
      cInf(4) = DCMPLX( 0.0d0,-aInf )
      cInf(5) = DCMPLX( aInf,  aInf )

      // .. NaN entries ..
      aNaN = aInf / aInf
      cNaN(1) = DCMPLX( aNaN, 0.0d0 )
      cNaN(2) = DCMPLX( 0.0d0, aNaN )
      cNaN(3) = DCMPLX( aNaN,  aNaN )


      // .. Tests ..

      // Test (a) Infs
      for (i = 1; i <= nInf; i++) { // 10
          nTests = nTests + 3
          Y = cInf(i)
          R = czero * Y
          if (R .eq. R) {
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'ia',i, czero, Y, R, 'NaN'
          }
          R = cone * Y
          if ( (R .ne. Y) .and. (R .eq. R) ) {
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'ib',i, cone, Y, R, 'the input and NaN'
          }
          R = Y * Y
          if ( (i.eq.1) .or. (i.eq.2) ) {
              if ( (R .ne. cInf(1)) .and. (R .eq. R) ) {
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, 'Inf and NaN'
              }
          } else if ( (i.eq.3) .or. (i.eq.4) ) {
              if ( (R .ne. cInf(2)) .and. (R .eq. R) ) {
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, '-Inf and NaN'
              }
          } else {
              if (R .eq. R) {
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, 'NaN'
              }
          }
      } // 10

      // Test (b) NaNs
      for (i = 1; i <= nNaN; i++) { // 20
          nTests = nTests + 3
          Y = cNaN(i)
          R = czero * Y
          if (R .eq. R) {
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'na',i, czero, Y, R, 'NaN'
          }
          R = cone * Y
          if (R .eq. R) {
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'nb',i, cone, Y, R, 'NaN'
          }
          R = Y * Y
          if (R .eq. R) {
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'nc',i, Y, Y, R, 'NaN'
          }
      } // 20

      if (nFailingTests .gt. 0) {
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests, " pass for complex multiplication,", nFailingTests," fail."
      } else {
         print *, "# All tests pass for complex multiplication."
      }

      // .. Formats ..
 9998 FORMAT( '[',A2,I1, '] (', (ES24.16E3,SP,ES24.16E3,"*I"), ') * (', (ES24.16E3,SP,ES24.16E3,"*I"), ') = (', (ES24.16E3,SP,ES24.16E3,"*I"), ') differs from ', A17 )

      // End of zmul

      }
