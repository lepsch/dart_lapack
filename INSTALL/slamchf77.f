      REAL FUNCTION SLAMCH( CMACH )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             CMACH;
      // ..
      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               FIRST, LRND;
      int                BETA, IMAX, IMIN, IT;
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, RND, SFMIN, SMALL, T
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAMC2
      // ..
      // .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, EMAX, RMAX, PREC
      // ..
      // .. Data statements ..
      DATA               FIRST / true /
      // ..
      // .. Executable Statements ..

      if ( FIRST ) {
         slamc2(BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX );
         BASE = BETA
         T = IT
         if ( LRND ) {
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         } else {
            RND = ZERO
            EPS = BASE**( 1-IT )
         }
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         if ( SMALL.GE.SFMIN ) {

            // Use SMALL plus a bit, to avoid the possibility of rounding
            // causing overflow when computing  1/sfmin.

            SFMIN = SMALL*( ONE+EPS )
         }
      }

      if ( LSAME( CMACH, 'E' ) ) {
         RMACH = EPS
      } else if ( LSAME( CMACH, 'S' ) ) {
         RMACH = SFMIN
      } else if ( LSAME( CMACH, 'B' ) ) {
         RMACH = BASE
      } else if ( LSAME( CMACH, 'P' ) ) {
         RMACH = PREC
      } else if ( LSAME( CMACH, 'N' ) ) {
         RMACH = T
      } else if ( LSAME( CMACH, 'R' ) ) {
         RMACH = RND
      } else if ( LSAME( CMACH, 'M' ) ) {
         RMACH = EMIN
      } else if ( LSAME( CMACH, 'U' ) ) {
         RMACH = RMIN
      } else if ( LSAME( CMACH, 'L' ) ) {
         RMACH = EMAX
      } else if ( LSAME( CMACH, 'O' ) ) {
         RMACH = RMAX
      }

      SLAMCH = RMACH
      FIRST  = false;
      RETURN

      // End of SLAMCH

      }

************************************************************************
*> \brief \b SLAMC1
*> \details
*> \b Purpose:
*> \verbatim
*> SLAMC1 determines the machine parameters given by BETA, T, RND, and
*> IEEE1.
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          The number of ( BETA ) digits in the mantissa.
*> \endverbatim
*>
*> \param[out] RND
*> \verbatim
*>          Specifies whether proper rounding  ( RND = true )  or
*>          chopping  ( RND = false )  occurs in addition. This may not
*>          be a reliable guide to the way in which the machine performs
*>          its arithmetic.
*> \endverbatim
*>
*> \param[out] IEEE1
*> \verbatim
*>          Specifies whether rounding appears to be done in the IEEE
*>          'round to nearest' style.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*>
*> \ingroup lamc1
*>
*> \details \b Further \b Details
*> \verbatim
*>
*>  The routine is based on the routine  ENVRON  by Malcolm and
*>  incorporates suggestions by Gentleman and Marovich. See
*>
*>     Malcolm M. A. (1972) Algorithms to reveal properties of
*>        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*>
*>     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*>        that reveal properties of floating point arithmetic units.
*>        Comms. of the ACM, 17, 276-277.
*> \endverbatim
*>
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      bool               IEEE1, RND;
      int                BETA, T;
      // ..
* =====================================================================

      // .. Local Scalars ..
      bool               FIRST, LIEEE1, LRND;
      int                LBETA, LT;
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
      // ..
      // .. External Functions ..
      REAL               SLAMC3
      // EXTERNAL SLAMC3
      // ..
      // .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
      // ..
      // .. Data statements ..
      DATA               FIRST / true /
      // ..
      // .. Executable Statements ..

      if ( FIRST ) {
         ONE = 1

         // LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
         // IEEE1, T and RND.

         // Throughout this routine  we use the function  SLAMC3  to ensure
         // that relevant values are  stored and not held in registers,  or
         // are not affected by optimizers.

         // Compute  a = 2.0**m  with the  smallest positive integer m such
         // that

            // fl( a + 1.0 ) = a.

         A = 1
         C = 1

*+       WHILE( C == ONE )LOOP
         } // 10
         if ( C == ONE ) {
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         }
*+       END WHILE

         // Now compute  b = 2.0**m  with the smallest positive integer m
         // such that

            // fl( a + b ) .gt. a.

         B = 1
         C = SLAMC3( A, B )

*+       WHILE( C == A )LOOP
         } // 20
         if ( C == A ) {
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         }
*+       END WHILE

         // Now compute the base.  a and c  are neighbouring floating point
         // numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
         // their difference is beta. Adding 0.25 to c is to ensure that it
         // is truncated to beta and not ( beta - 1 ).

         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR

         // Now determine whether rounding or chopping occurs,  by adding a
         // bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         if ( C == A ) {
            LRND = true;
         } else {
            LRND = false;
         }
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) && ( C == A ) ) LRND = false;

         // Try and decide whether rounding is done in the  IEEE  'round to
         // nearest' style. B/2 is half a unit in the last place of the two
         // numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
         // zero, and SAVEC is odd. Thus adding B/2 to A should not  change
         // A, but adding B/2 to SAVEC should change SAVEC.

         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1 == A ) && ( T2.GT.SAVEC ) && LRND

         // Now find  the  mantissa, t.  It should  be the  integer part of
         // log to the base beta of a,  however it is safer to determine  t
         // by powering.  So we find t as the smallest positive integer for
         // which

            // fl( beta**t + 1.0 ) = 1.0.

         LT = 0
         A = 1
         C = 1

*+       WHILE( C == ONE )LOOP
         } // 30
         if ( C == ONE ) {
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         }
*+       END WHILE

      }

      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = false;
      RETURN

      // End of SLAMC1

      }

************************************************************************
*> \brief \b SLAMC2
*> \details
*> \b Purpose:
*> \verbatim
*> SLAMC2 determines the machine parameters specified in its argument
*> list.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*>
*> \ingroup lamc2
*>
*> \param[out] BETA
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          The number of ( BETA ) digits in the mantissa.
*> \endverbatim
*>
*> \param[out] RND
*> \verbatim
*>          Specifies whether proper rounding  ( RND = true )  or
*>          chopping  ( RND = false )  occurs in addition. This may not
*>          be a reliable guide to the way in which the machine performs
*>          its arithmetic.
*> \endverbatim
*>
*> \param[out] EPS
*> \verbatim
*>          The smallest positive number such that
*>             fl( 1.0 - EPS ) .LT. 1.0,
*>          where fl denotes the computed value.
*> \endverbatim
*>
*> \param[out] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow occurs.
*> \endverbatim
*>
*> \param[out] RMIN
*> \verbatim
*>          The smallest normalized number for the machine, given by
*>          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*>          of BETA.
*> \endverbatim
*>
*> \param[out] EMAX
*> \verbatim
*>          The maximum exponent before overflow occurs.
*> \endverbatim
*>
*> \param[out] RMAX
*> \verbatim
*>          The largest positive number for the machine, given by
*>          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*>          value of BETA.
*> \endverbatim
*>
*> \details \b Further \b Details
*> \verbatim
*>
*>  The computation of  EPS  is based on a routine PARANOIA by
*>  W. Kahan of the University of California at Berkeley.
*> \endverbatim
*>
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      bool               RND;
      int                BETA, EMAX, EMIN, T;
      REAL               EPS, RMAX, RMIN
      // ..
* =====================================================================

      // .. Local Scalars ..
      bool               FIRST, IEEE, IWARN, LIEEE1, LRND;
      int                GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, NGNMIN, NGPMIN;
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, SIXTH, SMALL, THIRD, TWO, ZERO;
      // ..
      // .. External Functions ..
      REAL               SLAMC3
      // EXTERNAL SLAMC3
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAMC1, SLAMC4, SLAMC5
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, LRMIN, LT
      // ..
      // .. Data statements ..
      DATA               FIRST / true / , IWARN / false /
      // ..
      // .. Executable Statements ..

      if ( FIRST ) {
         ZERO = 0
         ONE = 1
         TWO = 2

         // LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
         // BETA, T, RND, EPS, EMIN and RMIN.

         // Throughout this routine  we use the function  SLAMC3  to ensure
         // that relevant values are stored  and not held in registers,  or
         // are not affected by optimizers.

         // SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

         slamc1(LBETA, LT, LRND, LIEEE1 );

         // Start to find EPS.

         B = LBETA
         A = B**( -LT )
         LEPS = A

         // Try some tricks to see whether or not this is the correct  EPS.

         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         if (B.LT.LEPS) B = LEPS;

         LEPS = 1

*+       WHILE( ( LEPS.GT.B ) && ( B.GT.ZERO ) )LOOP
         } // 10
         if ( ( LEPS.GT.B ) && ( B.GT.ZERO ) ) {
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         }
*+       END WHILE

         if (A.LT.LEPS) LEPS = A;

         // Computation of EPS complete.

         // Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
         // Keep dividing  A by BETA until (gradual) underflow occurs. This
         // is detected when we cannot recover the previous A.

         RBASE = ONE / LBETA
         SMALL = ONE
         for (I = 1; I <= 3; I++) { // 20
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
         } // 20
         A = SLAMC3( ONE, SMALL )
         slamc4(NGPMIN, ONE, LBETA );
         slamc4(NGNMIN, -ONE, LBETA );
         slamc4(GPMIN, A, LBETA );
         slamc4(GNMIN, -A, LBETA );
         IEEE = false;

         if ( ( NGPMIN == NGNMIN ) && ( GPMIN == GNMIN ) ) {
            if ( NGPMIN == GPMIN ) {
               LEMIN = NGPMIN
             // ( Non twos-complement machines, no gradual underflow;
               // e.g.,  VAX )
            } else if ( ( GPMIN-NGPMIN ) == 3 ) {
               LEMIN = NGPMIN - 1 + LT
               IEEE = true;
             // ( Non twos-complement machines, with gradual underflow;
               // e.g., IEEE standard followers )
            } else {
               LEMIN = MIN( NGPMIN, GPMIN )
             // ( A guess; no known machine )
               IWARN = true;
            }

         } else if ( ( NGPMIN == GPMIN ) && ( NGNMIN == GNMIN ) ) {
            if ( ABS( NGPMIN-NGNMIN ) == 1 ) {
               LEMIN = MAX( NGPMIN, NGNMIN )
             // ( Twos-complement machines, no gradual underflow;
               // e.g., CYBER 205 )
            } else {
               LEMIN = MIN( NGPMIN, NGNMIN )
             // ( A guess; no known machine )
               IWARN = true;
            }

         } else if ( ( ABS( NGPMIN-NGNMIN ) == 1 ) && ( GPMIN == GNMIN ) ) {
            if ( ( GPMIN-MIN( NGPMIN, NGNMIN ) ) == 3 ) {
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
             // ( Twos-complement machines with gradual underflow;
               // no known machine )
            } else {
               LEMIN = MIN( NGPMIN, NGNMIN )
             // ( A guess; no known machine )
               IWARN = true;
            }

         } else {
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
          // ( A guess; no known machine )
            IWARN = true;
         }
         FIRST = false;
***
* Comment out this if block if EMIN is ok
         if ( IWARN ) {
            FIRST = true;
            WRITE( 6, FMT = 9999 )LEMIN
         }
***

         // Assume IEEE arithmetic if we found denormalised  numbers above,
         // or if arithmetic seems to round in the  IEEE style,  determined
         // in routine SLAMC1. A true IEEE machine should have both  things
         // true; however, faulty machines may have one or the other.

         IEEE = IEEE || LIEEE1

         // Compute  RMIN by successive division by  BETA. We could compute
         // RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
         // this computation.

         LRMIN = 1
         for (I = 1; I <= 1 - LEMIN; I++) { // 30
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
         } // 30

         // Finally, call SLAMC5 to compute EMAX and RMAX.

         slamc5(LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX );
      }

      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX

      RETURN

 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', '  EMIN = ', I8, / ' If, after inspection, the value EMIN looks', ' acceptable please comment out ', / ' the IF block as marked within the code of routine', ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )

      // End of SLAMC2

      }

************************************************************************
*> \brief \b SLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> SLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*>
*> \param[in] A
*>
*> \param[in] B
*> \verbatim
*>          The values A and B.
*> \endverbatim
*>
*> \ingroup lamc3
*>
      REAL FUNCTION SLAMC3( A, B )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      REAL               A, B
      // ..
* =====================================================================

      // .. Executable Statements ..

      SLAMC3 = A + B

      RETURN

      // End of SLAMC3

      }

************************************************************************
*> \brief \b SLAMC4
*> \details
*> \b Purpose:
*> \verbatim
*> SLAMC4 is a service routine for SLAMC2.
*> \endverbatim
*>
*> \param[out] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow, computed by
*>          setting A = START and dividing by BASE until the previous A
*>          can not be recovered.
*> \endverbatim
*>
*> \param[in] START
*> \verbatim
*>          The starting point for determining EMIN.
*> \endverbatim
*>
*> \param[in] BASE
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
*> \ingroup lamc4
*>
      SUBROUTINE SLAMC4( EMIN, START, BASE )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      int                BASE;
      int                EMIN;
      REAL               START
      // ..
* =====================================================================

      // .. Local Scalars ..
      int                I;
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
      // ..
      // .. External Functions ..
      REAL               SLAMC3
      // EXTERNAL SLAMC3
      // ..
      // .. Executable Statements ..

      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1 == A ) && ( C2 == A ).AND.
*    $       ( D1 == A ) && ( D2 == A )      )LOOP
      } // 10
      if ( ( C1 == A ) && ( C2 == A ) && ( D1 == A ) && ( D2 == A ) ) {
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         for (I = 1; I <= BASE; I++) { // 20
            D1 = D1 + B1
         } // 20
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         for (I = 1; I <= BASE; I++) { // 30
            D2 = D2 + B2
         } // 30
         GO TO 10
      }
*+    END WHILE

      RETURN

      // End of SLAMC4

      }

************************************************************************
*> \brief \b SLAMC5
*> \details
*> \b Purpose:
*> \verbatim
*> SLAMC5 attempts to compute RMAX, the largest machine floating-point
*> number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*> approximately to a power of 2.  It will fail on machines where this
*> assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*> EMAX = 28718).  It will also fail if the value supplied for EMIN is
*> too large (i.e. too close to zero), probably with overflow.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          The base of floating-point arithmetic.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          The number of base BETA digits in the mantissa of a
*>          floating-point value.
*> \endverbatim
*>
*> \param[in] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow.
*> \endverbatim
*>
*> \param[in] IEEE
*> \verbatim
*>          A logical flag specifying whether or not the arithmetic
*>          system is thought to comply with the IEEE standard.
*> \endverbatim
*>
*> \param[out] EMAX
*> \verbatim
*>          The largest exponent before overflow
*> \endverbatim
*>
*> \param[out] RMAX
*> \verbatim
*>          The largest machine floating-point number.
*> \endverbatim
*>
*> \ingroup lamc5
*>
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      bool               IEEE;
      int                BETA, EMAX, EMIN, P;
      REAL               RMAX
      // ..
* =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP;
      REAL               OLDY, RECBAS, Y, Z
      // ..
      // .. External Functions ..
      REAL               SLAMC3
      // EXTERNAL SLAMC3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // First compute LEXP and UEXP, two powers of 2 that bound
      // abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
      // approximately to the bound that is closest to abs(EMIN).
      // (EMAX is the exponent of the required number RMAX).

      LEXP = 1
      EXBITS = 1
      } // 10
      TRY = LEXP*2
      if ( TRY.LE.( -EMIN ) ) {
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      }
      if ( LEXP == -EMIN ) {
         UEXP = LEXP
      } else {
         UEXP = TRY
         EXBITS = EXBITS + 1
      }

      // Now -LEXP is less than or equal to EMIN, and -UEXP is greater
      // than or equal to EMIN. EXBITS is the number of bits needed to
      // store the exponent.

      if ( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) {
         EXPSUM = 2*LEXP
      } else {
         EXPSUM = 2*UEXP
      }

      // EXPSUM is the exponent range, approximately equal to
      // EMAX - EMIN + 1 .

      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P

      // NBITS is the total number of bits needed to store a
      // floating-point number.

      if ( ( MOD( NBITS, 2 ) == 1 ) && ( BETA == 2 ) ) {

         // Either there are an odd number of bits used to store a
         // floating-point number, which is unlikely, or some bits are
         // not used in the representation of numbers, which is possible,
         // (e.g. Cray machines) or the mantissa has an implicit bit,
         // (e.g. IEEE machines, Dec Vax machines), which is perhaps the
         // most likely. We have to assume the last alternative.
         // If this is true, then we need to reduce EMAX by one because
         // there must be some way of representing zero in an implicit-bit
         // system. On machines like Cray, we are reducing EMAX by one
         // unnecessarily.

         EMAX = EMAX - 1
      }

      if ( IEEE ) {

         // Assume we are on an IEEE machine which reserves one exponent
         // for infinity and NaN.

         EMAX = EMAX - 1
      }

      // Now create RMAX, the largest machine number, which should
      // be equal to (1.0 - BETA**(-P)) * BETA**EMAX .

      // First compute 1.0 - BETA**(-P), being careful that the
      // result is less than 1.0 .

      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      for (I = 1; I <= P; I++) { // 20
         Z = Z*RECBAS
         if (Y.LT.ONE) OLDY = Y;
         Y = SLAMC3( Y, Z )
      } // 20
      if (Y.GE.ONE) Y = OLDY;

      // Now multiply by BETA**EMAX to get RMAX.

      for (I = 1; I <= EMAX; I++) { // 30
         Y = SLAMC3( Y*BETA, ZERO )
      } // 30

      RMAX = Y
      RETURN

      // End of SLAMC5

      }
