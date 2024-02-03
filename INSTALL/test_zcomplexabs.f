void main() {
*  -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // ..
      // .. Local parameters ..
      bool              debug;
      const           debug = false ;
      int               N, nNaN, nInf;
      const           N = 4, nNaN = 3, nInf = 5 ;
      double            threeFourth, fiveFourth, oneHalf;
      const           threeFourth = 3.0d0 / 4, fiveFourth = 5.0d0 / 4, oneHalf = 1.0d0 / 2 ;
      // ..
      // .. Local Variables ..
      int               i, min, Max, m, subnormalTreatedAs0, caseAFails, caseBFails, caseCFails, caseDFails, caseEFails, caseFFails, nFailingTests, nTests;
      double            X( N ), R, answerC, answerD, aInf, aNaN, relDiff, b, eps, blueMin, blueMax, Xj, stepX(N), limX(N);
      double complex    Y, cInf( nInf ), cNaN( nNaN )

      // .. Intrinsic Functions ..
      // intrinsic ABS, DBLE, RADIX, CEILING, TINY, DIGITS, SQRT, MAXEXPONENT, MINEXPONENT, FLOOR, HUGE, DCMPLX, EPSILON


      // .. Initialize error counts ..
      subnormalTreatedAs0 = 0
      caseAFails = 0
      caseBFails = 0
      caseCFails = 0
      caseDFails = 0
      caseEFails = 0
      caseFFails = 0
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

      // .. Vector X ..
      X(1) = TINY(0.0d0) * b**( DBLE(1-m) )
      X(2) = TINY(0.0d0)
      X(3) = HUGE(0.0d0)
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
      aInf = X(3) * 2
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
      if (Xj == 0.0d0) {
        subnormalTreatedAs0 = subnormalTreatedAs0 + 1
        if (debug || subnormalTreatedAs0 == 1) {
            print *, "!! fl( subnormal ) may be 0"
        }
      } else {
        for (i = 1; i <= N; i++) { // 100
            Xj = X(i)
            if (Xj == 0.0d0) {
                subnormalTreatedAs0 = subnormalTreatedAs0 + 1
                if (debug || subnormalTreatedAs0 == 1) {
                    print *, "!! fl( subnormal ) may be 0"
                }
            }
        } // 100
      }

      // Test (a) y = x + 0 * I, |y| = x
      for (i = 1; i <= N; i++) { // 10
        Xj = X(i)
        if (Xj == 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug || subnormalTreatedAs0 == 1) {
                print *, "!! [a] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj != limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( Xj, 0.0d0 )
                R = ABS( Y )
                if (R != Xj) {
                    caseAFails = caseAFails + 1
                    if (caseAFails == 1) {
                        print *, "!! Some ABS(x+0*I) differ from ABS(x)"
                    }
                    WRITE( 0, FMT = 9999 ) 'a',i, Xj, '(1+0*I)', R, Xj
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 10

      // Test (b) y = 0 + x * I, |y| = x
      for (i = 1; i <= N; i++) { // 20
        Xj = X(i)
        if (Xj == 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug || subnormalTreatedAs0 == 1) {
                print *, "!! [b] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj != limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( 0.0d0, Xj )
                R = ABS( Y )
                if (R != Xj) {
                    caseBFails = caseBFails + 1
                    if (caseBFails == 1) {
                        print *, "!! Some ABS(0+x*I) differ from ABS(x)"
                    }
                    WRITE( 0, FMT = 9999 ) 'b',i, Xj, '(0+1*I)', R, Xj
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 20

      // Test (c) y = (3/4)*x + x * I, |y| = (5/4)*x
      for (i = 1; i <= N; i++) { // 30
        if (i == 3) go to 30;
        if (i == 1) {
            Xj = 4*X(i)
        } else {
            Xj = X(i)
        }
        if (Xj == 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug || subnormalTreatedAs0 == 1) {
                print *, "!! [c] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj != limX(i) )
                nTests = nTests + 1
                answerC = fiveFourth * Xj
                Y = DCMPLX( threeFourth * Xj, Xj )
                R = ABS( Y )
                if (R != answerC) {
                    caseCFails = caseCFails + 1
                    if (caseCFails == 1) {
                        print *,  "!! Some ABS(x*(3/4+I)) differ from (5/4)*ABS(x)"
                    }
                    WRITE( 0, FMT = 9999 ) 'c',i, Xj, '(3/4+I)', R, answerC
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 30

      // Test (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2)
      for (i = 1; i <= N; i++) { // 40
        if (i == 1) {
            Xj = 2*X(i)
        } else {
            Xj = X(i)
        }
        if (Xj == 0.0d0) {
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if (debug || subnormalTreatedAs0 == 1) {
                print *, "!! [d] fl( subnormal ) may be 0"
            }
        } else {
            do while( Xj != limX(i) )
                answerD = (oneHalf * Xj) * SQRT(2.0d0)
                if (answerD == 0.0d0) {
                    subnormalTreatedAs0 = subnormalTreatedAs0 + 1
                    if (debug || subnormalTreatedAs0 == 1) {
                        print *, "!! [d] fl( subnormal ) may be 0"
                    }
                } else {
                    nTests = nTests + 1
                    Y = DCMPLX( oneHalf * Xj, oneHalf * Xj )
                    R = ABS( Y )
                    relDiff = ABS(R-answerD)/answerD
                    if ( relDiff .ge. (0.5*eps) ) {
                        caseDFails = caseDFails + 1
                        if (caseDFails == 1) {
                            print *,  "!! Some ABS(x*(1+I)) differ from sqrt(2)*ABS(x)"
                        }
                        WRITE( 0, FMT = 9999 ) 'd',i, (oneHalf*Xj), '(1+1*I)', R, answerD
                    }
                }
                Xj = Xj * stepX(i)
            }
        }
      } // 40

      // Test (e) Infs
      for (i = 1; i <= nInf; i++) { // 50
        nTests = nTests + 1
        Y = cInf(i)
        R = ABS( Y )
        if ( .not.(R .gt. HUGE(0.0d0)) ) {
            caseEFails = caseEFails + 1
            WRITE( *, FMT = 9997 ) 'i',i, Y, R
        }
      } // 50

      // Test (f) NaNs
      for (i = 1; i <= nNaN; i++) { // 60
        nTests = nTests + 1
        Y = cNaN(i)
        R = ABS( Y )
        if (R == R) {
            caseFFails = caseFFails + 1
            WRITE( *, FMT = 9998 ) 'n',i, Y, R
        }
      } // 60

      // If any test fails, displays a message
      nFailingTests = caseAFails + caseBFails + caseCFails + caseDFails + caseEFails + caseFFails
      if (nFailingTests .gt. 0) {
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests, " pass for ABS(a+b*I),", nFailingTests, " tests fail."
      } else {
         print *, "# All tests pass for ABS(a+b*I)"
      }

      // If anything was written to stderr, print the message
      if ( (caseAFails .gt. 0) || (caseBFails .gt. 0) || (caseCFails .gt. 0) || (caseDFails .gt. 0) ) {
         print *, "# Please check the failed ABS(a+b*I) in [stderr]"
      }

      // .. Formats ..
 9997 FORMAT( '[',A1,I1, '] ABS(', (ES8.1,SP,ES8.1,"*I"), ' ) = ', ES8.1, ' differs from Inf' )

 9998 FORMAT( '[',A1,I1, '] ABS(', (ES8.1,SP,ES8.1,"*I"), ' ) = ', ES8.1, ' differs from NaN' )

 9999 FORMAT( '[',A1,I1, '] ABS(', ES24.16E3, ' * ', A7, ' ) = ', ES24.16E3, ' differs from ', ES24.16E3 )

      // End of zabs

      }
