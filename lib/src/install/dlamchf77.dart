// ignore_for_file: type_annotate_public_apis

import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';

class _DlamchCache {
  var FIRST = true;
  var EPS = 0.0,
      SFMIN = 0.0,
      BASE = 0.0,
      T = 0.0,
      RND = 0.0,
      EMIN = 0.0,
      RMIN = 0.0,
      EMAX = 0.0,
      RMAX = 0.0,
      PREC = 0.0;
}

final _dlamchCache = _DlamchCache();

double dlamch(final String CMACH) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  final LRND = Box(false);
  final BETA = Box(0), IMAX = Box(0), IMIN = Box(0), IT = Box(0);
  double RMACH = 0, SMALL = 0;

  if (_dlamchCache.FIRST) {
    var EPS = Box(0.0), RMIN = Box(0.0), RMAX = Box(0.0);
    dlamc2(BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX);
    _dlamchCache.EPS = EPS.value;
    _dlamchCache.RMIN = RMIN.value;
    _dlamchCache.RMAX = RMAX.value;
    _dlamchCache.BASE = BETA.value.toDouble();
    _dlamchCache.T = IT.value.toDouble();
    if (LRND.value) {
      _dlamchCache.RND = ONE;
      _dlamchCache.EPS = pow(_dlamchCache.BASE, 1 - IT.value) / 2;
    } else {
      _dlamchCache.RND = ZERO;
      _dlamchCache.EPS = pow(_dlamchCache.BASE, 1 - IT.value).toDouble();
    }
    _dlamchCache.PREC = _dlamchCache.EPS * _dlamchCache.BASE;
    _dlamchCache.EMIN = IMIN.value.toDouble();
    _dlamchCache.EMAX = IMAX.value.toDouble();
    _dlamchCache.SFMIN = _dlamchCache.RMIN;
    SMALL = ONE / _dlamchCache.RMAX;
    if (SMALL >= _dlamchCache.SFMIN) {
      // Use SMALL plus a bit, to avoid the possibility of rounding
      // causing overflow when computing  1/sfmin.

      _dlamchCache.SFMIN = SMALL * (ONE + _dlamchCache.EPS);
    }
  }

  if (lsame(CMACH, 'E')) {
    RMACH = _dlamchCache.EPS;
  } else if (lsame(CMACH, 'S')) {
    RMACH = _dlamchCache.SFMIN;
  } else if (lsame(CMACH, 'B')) {
    RMACH = _dlamchCache.BASE;
  } else if (lsame(CMACH, 'P')) {
    RMACH = _dlamchCache.PREC;
  } else if (lsame(CMACH, 'N')) {
    RMACH = _dlamchCache.T;
  } else if (lsame(CMACH, 'R')) {
    RMACH = _dlamchCache.RND;
  } else if (lsame(CMACH, 'M')) {
    RMACH = _dlamchCache.EMIN;
  } else if (lsame(CMACH, 'U')) {
    RMACH = _dlamchCache.RMIN;
  } else if (lsame(CMACH, 'L')) {
    RMACH = _dlamchCache.EMAX;
  } else if (lsame(CMACH, 'O')) {
    RMACH = _dlamchCache.RMAX;
  }

  _dlamchCache.FIRST = false;
  return RMACH;
}

// ***********************************************************************
// > \brief \b DLAMC1
// > \details
// > \b Purpose:
// > \verbatim
// > DLAMC1 determines the machine parameters given by BETA, T, RND, and
// > IEEE1.
// > \endverbatim
// >
// > \param[out] BETA
// > \verbatim
// >          The base of the machine.
// > \endverbatim
// >
// > \param[out] T
// > \verbatim
// >          The number of ( BETA ) digits in the mantissa.
// > \endverbatim
// >
// > \param[out] RND
// > \verbatim
// >          Specifies whether proper rounding  ( RND = true )  or
// >          chopping  ( RND = false )  occurs in addition. This may not
// >          be a reliable guide to the way in which the machine performs
// >          its arithmetic.
// > \endverbatim
// >
// > \param[out] IEEE1
// > \verbatim
// >          Specifies whether rounding appears to be done in the IEEE
// >          'round to nearest' style.
// > \endverbatim
// > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
// >
// > \ingroup lamc1
// >
// > \details \b Further \b Details
// > \verbatim
// >
// >  The routine is based on the routine  ENVRON  by Malcolm and
// >  incorporates suggestions by Gentleman and Marovich. See
// >
// >     Malcolm M. A. (1972) Algorithms to reveal properties of
// >        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
// >
// >     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
// >        that reveal properties of floating point arithmetic units.
// >        Comms. of the ACM, 17, 276-277.
// > \endverbatim
// >

class _Dlamc1Cache {
  bool FIRST = true, LIEEE1 = false, LRND = false;
  int LBETA = 0, LT = 0;
}

final _dlamc1Cache = _Dlamc1Cache();

void dlamc1(
  final Box<int> BETA,
  final Box<int> T,
  final Box<bool> RND,
  final Box<bool> IEEE1,
) {
// -- LAPACK auxiliary routine --
// -- Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  double A, B, C, F, ONE, QTR, SAVEC, T1, T2;

  if (_dlamc1Cache.FIRST) {
    ONE = 1;

    // LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
    // IEEE1, T and RND.

    // Throughout this routine  we use the function  dlamc3  to ensure
    // that relevant values are  stored and not held in registers,  or
    // are not affected by optimizers.

    // Compute  a = 2.0**m  with the  smallest positive integer m such
    // that

    //    fl( a + 1.0 ) = a.

    A = 1;
    C = 1;

    while (C == ONE) {
      A = 2 * A;
      C = dlamc3(A, ONE);
      C = dlamc3(C, -A);
    }

    // Now compute  b = 2.0**m  with the smallest positive integer m
    // such that

    //    fl( a + b ) > a.

    B = 1;
    C = dlamc3(A, B);

    while (C == A) {
      B = 2 * B;
      C = dlamc3(A, B);
    }

    // Now compute the base.  a and c  are neighbouring floating point
    // numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
    // their difference is beta. Adding 0.25 to c is to ensure that it
    // is truncated to beta and not ( beta - 1 ).

    QTR = ONE / 4;
    SAVEC = C;
    C = dlamc3(C, -A);
    _dlamc1Cache.LBETA = (C + QTR).toInt();

    // Now determine whether rounding or chopping occurs,  by adding a
    // bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

    B = _dlamc1Cache.LBETA.toDouble();
    F = dlamc3(B / 2, -B / 100);
    C = dlamc3(F, A);
    if (C == A) {
      _dlamc1Cache.LRND = true;
    } else {
      _dlamc1Cache.LRND = false;
    }
    F = dlamc3(B / 2, B / 100);
    C = dlamc3(F, A);
    if ((_dlamc1Cache.LRND) && (C == A)) _dlamc1Cache.LRND = false;

    // Try and decide whether rounding is done in the  IEEE  'round to
    // nearest' style. B/2 is half a unit in the last place of the two
    // numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
    // zero, and SAVEC is odd. Thus adding B/2 to A should not  change
    // A, but adding B/2 to SAVEC should change SAVEC.

    T1 = dlamc3(B / 2, A);
    T2 = dlamc3(B / 2, SAVEC);
    _dlamc1Cache.LIEEE1 = (T1 == A) && (T2 > SAVEC) && _dlamc1Cache.LRND;

    // Now find  the  mantissa, t.  It should  be the  integer part of
    // log to the base beta of a,  however it is safer to determine  t
    // by powering.  So we find t as the smallest positive integer for
    // which

    // fl( beta**t + 1.0 ) = 1.0.

    _dlamc1Cache.LT = 0;
    A = 1;
    C = 1;

    while (C == ONE) {
      _dlamc1Cache.LT += 1;
      A = A * _dlamc1Cache.LBETA;
      C = dlamc3(A, ONE);
      C = dlamc3(C, -A);
    }
  }

  BETA.value = _dlamc1Cache.LBETA;
  T.value = _dlamc1Cache.LT;
  RND.value = _dlamc1Cache.LRND;
  IEEE1.value = _dlamc1Cache.LIEEE1;
  _dlamc1Cache.FIRST = false;
}

// ***********************************************************************
// > \brief \b DLAMC2
// > \details
// > \b Purpose:
// > \verbatim
// > DLAMC2 determines the machine parameters specified in its argument
// > list.
// > \endverbatim
// > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
// >
// > \ingroup lamc2
// >
// > \param[out] BETA
// > \verbatim
// >          The base of the machine.
// > \endverbatim
// >
// > \param[out] T
// > \verbatim
// >          The number of ( BETA ) digits in the mantissa.
// > \endverbatim
// >
// > \param[out] RND
// > \verbatim
// >          Specifies whether proper rounding  ( RND = true )  or
// >          chopping  ( RND = false )  occurs in addition. This may not
// >          be a reliable guide to the way in which the machine performs
// >          its arithmetic.
// > \endverbatim
// >
// > \param[out] EPS
// > \verbatim
// >          The smallest positive number such that
// >             fl( 1.0 - EPS ) < 1.0,
// >          where fl denotes the computed value.
// > \endverbatim
// >
// > \param[out] EMIN
// > \verbatim
// >          The minimum exponent before (gradual) underflow occurs.
// > \endverbatim
// >
// > \param[out] RMIN
// > \verbatim
// >          The smallest normalized number for the machine, given by
// >          BASE**( EMIN - 1 ), where  BASE  is the floating point value
// >          of BETA.
// > \endverbatim
// >
// > \param[out] EMAX
// > \verbatim
// >          The maximum exponent before overflow occurs.
// > \endverbatim
// >
// > \param[out] RMAX
// > \verbatim
// >          The largest positive number for the machine, given by
// >          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
// >          value of BETA.
// > \endverbatim
// >
// > \details \b Further \b Details
// > \verbatim
// >
// >  The computation of  EPS  is based on a routine PARANOIA by
// >  W. Kahan of the University of California at Berkeley.
// > \endverbatim

class _Dlamc2Cache {
  bool FIRST = true, IWARN = false;
  int LBETA = 0, LEMAX = 0, LEMIN = 0, LT = 0;
  double LEPS = 0.0, LRMAX = 0.0, LRMIN = 0.0;
}

final _dlamc2Cache = _Dlamc2Cache();

void dlamc2(
  final Box<int> BETA,
  final Box<int> T,
  final Box<bool> RND,
  final Box<double> EPS,
  final Box<int> EMIN,
  final Box<double> RMIN,
  final Box<int> EMAX,
  final Box<double> RMAX,
) {
// -- LAPACK auxiliary routine --
// -- Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  bool IEEE, IWARN = false;
  int I, LEMIN = 0;
  double A,
      B,
      C,
      HALF,
      LEPS = 0,
      LRMIN = 0,
      ONE,
      RBASE,
      SIXTH,
      SMALL,
      THIRD,
      TWO,
      ZERO;
  final LRND = Box(false), LIEEE1 = Box(false);
  final LEMAX = Box(0),
      NGPMIN = Box(0),
      NGNMIN = Box(0),
      GNMIN = Box(0),
      GPMIN = Box(0),
      LBETA = Box(0),
      LT = Box(0);
  final LRMAX = Box(0.0);

  if (_dlamc2Cache.FIRST) {
    ZERO = 0;
    ONE = 1;
    TWO = 2;

    // LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
    // BETA, T, RND, EPS, EMIN and RMIN.

    // Throughout this routine  we use the function  dlamc3  to ensure
    // that relevant values are stored  and not held in registers,  or
    // are not affected by optimizers.

    // DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

    dlamc1(LBETA, LT, LRND, LIEEE1);

    // Start to find EPS.

    B = LBETA.value.toDouble();
    A = pow(B, -LT.value).toDouble();
    LEPS = A;

    // Try some tricks to see whether or not this is the correct  EPS.

    B = TWO / 3;
    HALF = ONE / 2;
    SIXTH = dlamc3(B, -HALF);
    THIRD = dlamc3(SIXTH, SIXTH);
    B = dlamc3(THIRD, -HALF);
    B = dlamc3(B, SIXTH);
    B = (B).abs();
    if (B < LEPS) B = LEPS;

    LEPS = 1;

    while ((LEPS > B) && (B > ZERO)) {
      LEPS = B;
      C = dlamc3(HALF * LEPS, pow(TWO, 5) * pow(LEPS, 2).toDouble());
      C = dlamc3(HALF, -C);
      B = dlamc3(HALF, C);
      C = dlamc3(HALF, -B);
      B = dlamc3(HALF, C);
    }

    if (A < LEPS) LEPS = A;

    // Computation of EPS complete.

    // Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
    // Keep dividing  A by BETA until (gradual) underflow occurs. This
    // is detected when we cannot recover the previous A.

    RBASE = ONE / LBETA.value;
    SMALL = ONE;
    for (I = 1; I <= 3; I++) {
      SMALL = dlamc3(SMALL * RBASE, ZERO);
    }
    A = dlamc3(ONE, SMALL);
    dlamc4(NGPMIN, ONE, LBETA.value);
    dlamc4(NGNMIN, -ONE, LBETA.value);
    dlamc4(GPMIN, A, LBETA.value);
    dlamc4(GNMIN, -A, LBETA.value);
    IEEE = false;

    if ((NGPMIN == NGNMIN) && (GPMIN == GNMIN)) {
      if (NGPMIN == GPMIN) {
        LEMIN = NGPMIN.value;
        // ( Non twos-complement machines, no gradual underflow;
        // e.g.,  VAX )
      } else if ((GPMIN.value - NGPMIN.value) == 3) {
        LEMIN = NGPMIN.value - 1 + LT.value;
        IEEE = true;
        // ( Non twos-complement machines, with gradual underflow;
        // e.g., IEEE standard followers )
      } else {
        LEMIN = min(NGPMIN.value, GPMIN.value);
        // ( A guess; no known machine )
        IWARN = true;
      }
    } else if ((NGPMIN == GPMIN) && (NGNMIN == GNMIN)) {
      if ((NGPMIN.value - NGNMIN.value).abs() == 1) {
        LEMIN = max(NGPMIN.value, NGNMIN.value);
        // ( Twos-complement machines, no gradual underflow;
        // e.g., CYBER 205 )
      } else {
        LEMIN = min(NGPMIN.value, NGNMIN.value);
        // ( A guess; no known machine )
        IWARN = true;
      }
    } else if (((NGPMIN.value - NGNMIN.value).abs() == 1) && (GPMIN == GNMIN)) {
      if ((GPMIN.value - min(NGPMIN.value, NGNMIN.value)) == 3) {
        LEMIN = max(NGPMIN.value, NGNMIN.value) - 1 + LT.value;
        // ( Twos-complement machines with gradual underflow;
        // no known machine )
      } else {
        LEMIN = min(NGPMIN.value, NGNMIN.value);
        // ( A guess; no known machine )
        IWARN = true;
      }
    } else {
      LEMIN =
          min(min(NGPMIN.value, NGNMIN.value), min(GPMIN.value, GNMIN.value));
      // ( A guess; no known machine )
      IWARN = true;
    }
    _dlamc2Cache.FIRST = false;
    // **
    // Comment out this if block if EMIN is ok
    if (IWARN) {
      _dlamc2Cache.FIRST = true;
      print(
          '\n\n WARNING. The value EMIN may be incorrect:-  EMIN = ${LEMIN.i8}\n If, after inspection, the value EMIN looks acceptable please comment out \n the IF block as marked within the code of routine DLAMC2,\n otherwise supply EMIN explicitly.\n');
    }
    // **

    // Assume IEEE arithmetic if we found denormalised  numbers above,
    // or if arithmetic seems to round in the  IEEE style,  determined
    // in routine DLAMC1. A true IEEE machine should have both  things
    // true; however, faulty machines may have one or the other.

    IEEE = IEEE || LIEEE1.value;

    // Compute  RMIN by successive division by  BETA. We could compute
    // RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
    // this computation.

    LRMIN = 1;
    for (I = 1; I <= 1 - LEMIN; I++) {
      LRMIN = dlamc3(LRMIN * RBASE, ZERO);
    }

    // Finally, call DLAMC5 to compute EMAX and RMAX.

    dlamc5(LBETA.value, LT.value, LEMIN, IEEE, LEMAX, LRMAX);
  }

  BETA.value = LBETA.value;
  T.value = LT.value;
  RND.value = LRND.value;
  EPS.value = LEPS;
  EMIN.value = LEMIN;
  RMIN.value = LRMIN;
  EMAX.value = LEMAX.value;
  RMAX.value = LRMAX.value;
}

// ***********************************************************************
// > \brief \b dlamc3
// > \details
// > \b Purpose:
// > \verbatim
// > dlamc3  is intended to force  A  and  B  to be stored prior to doing
// > the addition of  A  and  B ,  for use in situations where optimizers
// > might hold one of these in a register.
// > \endverbatim
// >
// > \param[in] A
// >
// > \param[in] B
// > \verbatim
// >          The values A and B.
// > \endverbatim
// >
// > \ingroup lamc3
// >
double dlamc3(final double A, final double B) {
// -- LAPACK auxiliary routine --
// -- Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  return A + B;
}

// ***********************************************************************
// > \brief \b DLAMC4
// > \details
// > \b Purpose:
// > \verbatim
// > DLAMC4 is a service routine for DLAMC2.
// > \endverbatim
// >
// > \param[out] EMIN
// > \verbatim
// >          The minimum exponent before (gradual) underflow, computed by
// >          setting A = START and dividing by BASE until the previous A
// >          can not be recovered.
// > \endverbatim
// >
// > \param[in] START
// > \verbatim
// >          The starting point for determining EMIN.
// > \endverbatim
// >
// > \param[in] BASE
// > \verbatim
// >          The base of the machine.
// > \endverbatim
// >
// > \ingroup lamc4
// >
void dlamc4(final Box<int> EMIN, final double START, final int BASE) {
// -- LAPACK auxiliary routine --
// -- Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  int I;
  double A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO;

  A = START;
  ONE = 1;
  RBASE = ONE / BASE;
  ZERO = 0;
  EMIN.value = 1;
  B1 = dlamc3(A * RBASE, ZERO);
  C1 = A;
  C2 = A;
  D1 = A;
  D2 = A;
  while ((C1 == A) && (C2 == A) && (D1 == A) && (D2 == A)) {
    EMIN.value -= 1;
    A = B1;
    B1 = dlamc3(A / BASE, ZERO);
    C1 = dlamc3(B1 * BASE, ZERO);
    D1 = ZERO;
    for (I = 1; I <= BASE; I++) {
      D1 += B1;
    }
    B2 = dlamc3(A * RBASE, ZERO);
    C2 = dlamc3(B2 / RBASE, ZERO);
    D2 = ZERO;
    for (I = 1; I <= BASE; I++) {
      D2 += B2;
    }
  }
}

// ***********************************************************************
// > \brief \b DLAMC5
// > \details
// > \b Purpose:
// > \verbatim
// > DLAMC5 attempts to compute RMAX, the largest machine floating-point
// > number, without overflow.  It assumes that EMAX + abs(EMIN) sum
// > approximately to a power of 2.  It will fail on machines where this
// > assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
// > EMAX = 28718).  It will also fail if the value supplied for EMIN is
// > too large (i.e. too close to zero), probably with overflow.
// > \endverbatim
// >
// > \param[in] BETA
// > \verbatim
// >          The base of floating-point arithmetic.
// > \endverbatim
// >
// > \param[in] P
// > \verbatim
// >          The number of base BETA digits in the mantissa of a
// >          floating-point value.
// > \endverbatim
// >
// > \param[in] EMIN
// > \verbatim
// >          The minimum exponent before (gradual) underflow.
// > \endverbatim
// >
// > \param[in] IEEE
// > \verbatim
// >          A logical flag specifying whether or not the arithmetic
// >          system is thought to comply with the IEEE standard.
// > \endverbatim
// >
// > \param[out] EMAX
// > \verbatim
// >          The largest exponent before overflow
// > \endverbatim
// >
// > \param[out] RMAX
// > \verbatim
// >          The largest machine floating-point number.
// > \endverbatim
// >
// > \ingroup lamc5
// >
void dlamc5(
  final int BETA,
  final int P,
  final int EMIN,
  final bool IEEE,
  final Box<int> EMAX,
  final Box<double> RMAX,
) {
// -- LAPACK auxiliary routine --
// -- Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  const ZERO = 0.0, ONE = 1.0;
  int EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP;
  double OLDY = 0, RECBAS, Y, Z;

  // First compute LEXP and UEXP, two powers of 2 that bound
  // abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
  // approximately to the bound that is closest to abs(EMIN).
  // (EMAX is the exponent of the required number RMAX).

  LEXP = 1;
  EXBITS = 1;
  TRY = LEXP * 2;
  while (TRY <= -EMIN) {
    LEXP = TRY;
    EXBITS++;
    TRY = LEXP * 2;
  }
  if (LEXP == -EMIN) {
    UEXP = LEXP;
  } else {
    UEXP = TRY;
    EXBITS++;
  }

  // Now -LEXP is less than or equal to EMIN, and -UEXP is greater
  // than or equal to EMIN. EXBITS is the number of bits needed to
  // store the exponent.

  if ((UEXP + EMIN) > (-LEXP - EMIN)) {
    EXPSUM = 2 * LEXP;
  } else {
    EXPSUM = 2 * UEXP;
  }

  // EXPSUM is the exponent range, approximately equal to
  // EMAX - EMIN + 1 .

  EMAX.value = EXPSUM + EMIN - 1;
  NBITS = 1 + EXBITS + P;

  // NBITS is the total number of bits needed to store a
  // floating-point number.

  if (((NBITS % 2) == 1) && (BETA == 2)) {
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

    EMAX.value -= 1;
  }

  if (IEEE) {
    // Assume we are on an IEEE machine which reserves one exponent
    // for infinity and NaN.

    EMAX.value -= 1;
  }

  // Now create RMAX, the largest machine number, which should
  // be equal to (1.0 - BETA**(-P)) * BETA**EMAX .

  // First compute 1.0 - BETA**(-P), being careful that the
  // result is less than 1.0 .

  RECBAS = ONE / BETA;
  Z = BETA - ONE;
  Y = ZERO;
  for (I = 1; I <= P; I++) {
    Z = Z * RECBAS;
    if (Y < ONE) OLDY = Y;
    Y = dlamc3(Y, Z);
  }
  if (Y >= ONE) Y = OLDY;

  // Now multiply by BETA**EMAX to get RMAX.

  for (I = 1; I <= EMAX.value; I++) {
    Y = dlamc3(Y * BETA, ZERO);
  }

  RMAX.value = Y;
}
