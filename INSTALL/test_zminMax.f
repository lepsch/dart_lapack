      program zmul

*  -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // ..
      // .. Parameters ..
      int               n;
      const           n = 8 ;
      double            zero;
      const           zero = 0.0d0 ;
      // ..
      // .. Local Variables ..
      int               i, nFailingTests, nTests;
      double            aInf, aNaN, OV, R, X(n), Y(n);

      // .. Intrinsic Functions ..
      // intrinsic HUGE, MIN, MAX


      // .. Initialize error counts ..
      nFailingTests = 0
      nTests = 0

      // .. Inf and NaN entries ..
      OV = HUGE(0.0d0)
      aInf = OV * 2
      aNaN = aInf / aInf
      X = (/ -aInf, zero, -aInf,  zero, aInf,  aInf, zero, aNaN /)
      Y = (/  zero, aInf,  aInf, -aInf, zero, -aInf, aNaN, zero /)


      // .. Tests ..

      do 10 i = 1, 3
          nTests = nTests + 2
          R = MIN( X(i), Y(i) )
          if( R .ne. X(i) ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R
          endif
          R = MAX( X(i), Y(i) )
          if( R .ne. Y(i) ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R
          endif
  10  continue
      do 20 i = 4, 6
          nTests = nTests + 2
          R = MIN( X(i), Y(i) )
          if( R .ne. Y(i) ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R
          endif
          R = MAX( X(i), Y(i) )
          if( R .ne. X(i) ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R
          endif
  20  continue
      do 30 i = 7, 8
          nTests = nTests + 2
          R = MIN( X(i), Y(i) )
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MIN', X(i), Y(i), R
          endif
          R = MAX( X(i), Y(i) )
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'i',i, 'MAX', X(i), Y(i), R
          endif
  30  continue

      if( nFailingTests .gt. 0 ) then
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests, " pass for intrinsic MIN and MAX,", nFailingTests," fail."
      else
         print *, "# All tests pass for intrinsic MIN and MAX."
      endif

      // .. Formats ..
 9998 FORMAT( '[',A1,I1, '] ', A3, '(', F5.0, ',', F5.0, ') = ', F5.0 )

      // End of zmul

      }
