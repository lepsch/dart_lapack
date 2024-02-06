      void main() {
// -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // ..
      // .. Parameters ..
      int               n;
      const           n = 8 ;
      double            zero;
      const           zero = 0.0 ;
      // ..
      // .. Local Variables ..
      int               i, nFailingTests, nTests;
      double            aInf, aNaN, OV, R, X(n), Y(n);

      // .. Intrinsic Functions ..
      // intrinsic HUGE, MIN, MAX


      // .. Initialize error counts ..
      nFailingTests = 0;
      nTests = 0;

      // .. Inf and NaN entries ..
      OV = HUGE(0.0);
      aInf = OV * 2;
      aNaN = aInf / aInf;
      X = (/ -aInf, zero, -aInf,  zero, aInf,  aInf, zero, aNaN /);
      Y = (/  zero, aInf,  aInf, -aInf, zero, -aInf, aNaN, zero /);


      // .. Tests ..

      for (i = 1; i <= 3; i++) { // 10
          nTests = nTests + 2;
          R = min( X(i), Y(i) );
          if ( R != X(i) ) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R;
          }
          R = max( X(i), Y(i) );
          if ( R != Y(i) ) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R;
          }
      } // 10
      for (i = 4; i <= 6; i++) { // 20
          nTests = nTests + 2;
          R = min( X(i), Y(i) );
          if ( R != Y(i) ) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R;
          }
          R = max( X(i), Y(i) );
          if ( R != X(i) ) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R;
          }
      } // 20
      for (i = 7; i <= 8; i++) { // 30
          nTests = nTests + 2;
          R = min( X(i), Y(i) );
          if (R == R) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R;
          }
          R = max( X(i), Y(i) );
          if (R == R) {
              nFailingTests = nFailingTests + 1;
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R;
          }
      } // 30

      if (nFailingTests > 0) {
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests, " pass for intrinsic MIN and MAX,", nFailingTests," fail.";
      } else {
         print *, "# All tests pass for intrinsic MIN and MAX.";
      }

      // .. Formats ..
 9998 FORMAT( '[',A1,I1, '] ${.a3}(${.f5_0},${.f5_0}) = ${.f5_0}');
      }
