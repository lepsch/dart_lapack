void main() {
*  -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // ..
      // .. Local parameters ..
      bool              debug;
      const           debug = .false. ;
      int               N, nNaN, nInf;
      const           N = 4, nNaN = 3, nInf = 5 ;
      double            threeFourth, fiveFourth;
      const           threeFourth = 3.0d0 / 4, fiveFourth = 5.0d0 / 4 ;
      double complex    czero, cone
      const           czero = DCMPLX( 0.0d0, 0.0d0 ), cone  = DCMPLX( 1.0d0, 0.0d0 ) ;
      // ..
      // .. Local Variables ..
      int               i, min, Max, m, subnormalTreatedAs0, caseAFails, caseBFails, caseCFails, caseDFails, caseEFails, caseFFails, caseInfFails, caseNaNFails, nFailingTests, nTests;
      double            X( N ), aInf, aNaN, b, eps, blueMin, blueMax, OV, Xj, stepX(N), limX(N);
      double complex    Y, Y2, R, cInf( nInf ), cNaN( nNaN )

      // .. Intrinsic Functions ..
      // intrinsic DCONJG, DBLE, RADIX, CEILING, TINY, DIGITS, MAXEXPONENT, MINEXPONENT, FLOOR, HUGE, DCMPLX, EPSILON


      // .. Initialize error counts ..
      subnormalTreatedAs0 = 0
      caseAFails = 0
      caseBFails = 0
      caseCFails = 0
      caseDFails = 0
      caseEFails = 0
      caseFFails = 0
      caseInfFails = 0
      caseNaNFails = 0
      nFailingTests = 0
      nTests = 0

      // .. Initialize machine constants ..
      min = MINEXPONENT(0.0d0)
      Max = MAXEXPONENT(0.0d0)
      m = DIGITS(0.0d0)
      b = DBLE(RADIX(0.0d0))
      eps = EPSILON(0.0d0)
      blueMin = b**CEILING( (min - 1) * 0.5d0 )
      blueMax = b**FLOOR( (Max - m + 1) * 0.5d0 )
      OV = HUGE(0.0d0)

      // .. Vector X ..
      X(1) = TINY(0.0d0) * b**( DBLE(1-m) )
      X(2) = TINY(0.0d0)
      X(3) = OV
      X(4) = b**( DBLE(Max-1) )

      // .. Then modify X using the step ..
      stepX(1) = 2.0
      stepX(2) = 2.0
      stepX(3) = 0.0
      stepX(4) = 0.5

      // .. Up to the value ..
      limX(1) = X(2)
      limX(2) = 1.0
      limX(3) = 0.0
      limX(4) = 2.0

      // .. Inf entries ..
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

      if (debug) {
        print *, '# X :=', X
        print *, '# Blue min constant :=', blueMin
        print *, '# Blue max constant :=', blueMax
      }

      Xj = X(1)
      if (Xj .eq. 0.0d0) {
        subnormalTreatedAs0 = subnormalTreatedAs0 + 1
        if (debug .or. subnormalTreatedAs0 .eq. 1) {
            print *, "!! fl( subnormal ) may be 0"
        }
      } else {
        for (i = 1; i <= N; i++) { // 100
            Xj = X(i)
            if (Xj .eq. 0.0d0) {
                subnormalTreatedAs0 = subnormalTreatedAs0 + 1
                if (debug .or. subnormalTreatedAs0 .eq. 1) {
                    print *, "!! fl( subnormal ) may be 0"
                }
            }
        } // 100
      }

      // Test (a) y = x + 0 * I, y/y = 1
      for (i = 1; i <= N; i++) { // 10
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [a] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( Xj, 0.0d0 )
                R = Y / Y
                if (R .ne. 1.0D0) {
                    caseAFails = caseAFails + 1
                    if (caseAFails .eq. 1) {
                        print *, "!! Some (x+0*I)/(x+0*I) differ from 1"
                    }
                    WRITE( 0, FMT = 9999 ) 'a',i, Xj, '(x+0*I)/(x+0*I)', R, 1.0D0
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 10

      // Test (b) y = 0 + x * I, y/y = 1
      for (i = 1; i <= N; i++) { // 20
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [b] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( 0.0d0, Xj )
                R = Y / Y
                if (R .ne. 1.0D0) {
                    caseBFails = caseBFails + 1
                    if (caseBFails .eq. 1) {
                        print *, "!! Some (0+x*I)/(0+x*I) differ from 1"
                    }
                    WRITE( 0, FMT = 9999 ) 'b',i, Xj, '(0+x*I)/(0+x*I)', R, 1.0D0
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 20

      // Test (c) y = x + x * I, y/y = 1
      for (i = 1; i <= N; i++) { // 30
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [c] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( Xj, Xj )
                R = Y / Y
                if (R .ne. 1.0D0) {
                    caseCFails = caseCFails + 1
                    if (caseCFails .eq. 1) {
                        print *, "!! Some (x+x*I)/(x+x*I) differ from 1"
                    }
                    WRITE( 0, FMT = 9999 ) 'c',i, Xj, '(x+x*I)/(x+x*I)', R, 1.0D0
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 30

      // Test (d) y1 = 0 + x * I, y2 = x + 0 * I, y1/y2 = I
      for (i = 1; i <= N; i++) { // 40
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [d] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y  = DCMPLX( 0.0d0, Xj )
                Y2 = DCMPLX( Xj, 0.0d0 )
                R = Y / Y2
                if ( R .ne. DCMPLX(0.0D0,1.0D0) ) {
                    caseDFails = caseDFails + 1
                    if (caseDFails .eq. 1) {
                        print *, "!! Some (0+x*I)/(x+0*I) differ from I"
                    }
                    WRITE( 0, FMT = 9999 ) 'd',i, Xj, '(0+x*I)/(x+0*I)', R, DCMPLX(0.0D0,1.0D0)
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 40

      // Test (e) y1 = 0 + x * I, y2 = x + 0 * I, y2/y1 = -I
      for (i = 1; i <= N; i++) { // 50
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [e] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y  = DCMPLX( 0.0d0, Xj )
                Y2 = DCMPLX( Xj, 0.0d0 )
                R = Y2 / Y
                if ( R .ne. DCMPLX(0.0D0,-1.0D0) ) {
                    caseEFails = caseEFails + 1
                    if (caseEFails .eq. 1) {
                        print *,"!! Some (x+0*I)/(0+x*I) differ from -I"
                    }
                    WRITE( 0, FMT = 9999 ) 'e',i, Xj, '(x+0*I)/(0+x*I)', R, DCMPLX(0.0D0,-1.0D0)
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 50

      // Test (f) y = x + x * I, y/conj(y) = I
      for (i = 1; i <= N; i++) { // 60
        Xj = X(i)
        if (Xj .eq. 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug .or. subnormalTreatedAs0 .eq. 1) {
                print *, "!! [f] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y  = DCMPLX( Xj, Xj )
                R = Y / DCONJG( Y )
                if ( R .ne. DCMPLX(0.0D0,1.0D0) ) {
                    caseFFails = caseFFails + 1
                    if (caseFFails .eq. 1) {
                        print *, "!! Some (x+x*I)/(x-x*I) differ from I"
                    }
                    WRITE( 0, FMT = 9999 ) 'f',i, Xj, '(x+x*I)/(x-x*I)', R, DCMPLX(0.0D0,1.0D0)
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 60

      // Test (g) Infs
      for (i = 1; i <= nInf; i++) { // 70
          nTests = nTests + 3
          Y = cInf(i)
          R = czero / Y
          if ( (R .ne. czero) .and. (R .eq. R) ) {
              caseInfFails = caseInfFails + 1
              WRITE( *, FMT = 9998 ) 'ia',i, czero, Y, R, 'NaN and 0'
          }
          R = cone / Y
          if ( (R .ne. czero) .and. (R .eq. R) ) {
              caseInfFails = caseInfFails + 1
              WRITE( *, FMT = 9998 ) 'ib',i, cone, Y, R, 'NaN and 0'
          }
          R = Y / Y
          if (R .eq. R) {
              caseInfFails = caseInfFails + 1
              WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, 'NaN'
          }
      } // 70

      // Test (h) NaNs
      for (i = 1; i <= nNaN; i++) { // 80
          nTests = nTests + 3
          Y = cNaN(i)
          R = czero / Y
          if (R .eq. R) {
              caseNaNFails = caseNaNFails + 1
              WRITE( *, FMT = 9998 ) 'na',i, czero, Y, R, 'NaN'
          }
          R = cone / Y
          if (R .eq. R) {
              caseNaNFails = caseNaNFails + 1
              WRITE( *, FMT = 9998 ) 'nb',i, cone, Y, R, 'NaN'
          }
          R = Y / Y
          if (R .eq. R) {
              caseNaNFails = caseNaNFails + 1
              WRITE( *, FMT = 9998 ) 'nc',i, Y, Y, R, 'NaN'
          }
      } // 80

      // If any test fails, displays a message
      nFailingTests = caseAFails + caseBFails + caseCFails + caseDFails + caseEFails + caseFFails + caseInfFails + caseNaNFails
      if (nFailingTests .gt. 0) {
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests, " pass for complex division,", nFailingTests," fail."
      } else {
         print *, "# All tests pass for complex division."
      }

      // If anything was written to stderr, print the message
      if( (caseAFails .gt. 0) .or. (caseBFails .gt. 0) .or. (caseCFails .gt. 0) .or. (caseDFails .gt. 0) .or. (caseEFails .gt. 0) .or. (caseFFails .gt. 0) ) print *, "# Please check the failed divisions in [stderr]"

      // .. Formats ..
 9998 FORMAT( '[',A2,I1, '] ', (ES24.16E3,SP,ES24.16E3,"*I"), ' * ', (ES24.16E3,SP,ES24.16E3,"*I"), ' = ', (ES24.16E3,SP,ES24.16E3,"*I"), ' differs from ', A10 )

 9999 FORMAT( '[',A2,I1, '] X = ', ES24.16E3, ' : ', A15, ' = ', (ES24.16E3,SP,ES24.16E3,"*I"), ' differs from ', (ES24.16E3,SP,ES24.16E3,"*I") )

      // End of zdiv

      }
