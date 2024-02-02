      SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0E+0, ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               ALPHI, ALPHR, BETA, BIGNUM, EPS, SMLNUM, XNORM
      COMPLEX            SAVEALPHA
*     ..
*     .. External Functions ..
      REAL               SCNRM2, SLAMCH, SLAPY3, SLAPY2
      COMPLEX            CLADIV
      EXTERNAL           SCNRM2, SLAMCH, SLAPY3, SLAPY2, CLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSCAL, CSSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      EPS = SLAMCH( 'Precision' )
      XNORM = SCNRM2( N-1, X, INCX )
      ALPHR = REAL( ALPHA )
      ALPHI = AIMAG( ALPHA )
*
      IF( XNORM.LE.EPS*ABS(ALPHA) .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.
*
         IF( ALPHR.GE.ZERO ) THEN
*           When TAU.eq.ZERO, the vector is special-cased to be
*           all zeros in the application routines.  We do not need
*           to clear it.
            TAU = ZERO
         ELSE
*           However, the application routines rely on explicit
*           zero checks when TAU.ne.ZERO, and we must clear X.
            TAU = TWO
            DO J = 1, N-1
               X( 1 + (J-1)*INCX ) = ZERO
            END DO
            ALPHA = -ALPHA
         END IF
      ELSE
*
*        general case
*
         BETA = SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'E' )
         BIGNUM = ONE / SMLNUM
*
         KNT = 0
         IF( ABS( BETA ).LT.SMLNUM ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
   10       CONTINUE
            KNT = KNT + 1
            CALL CSSCAL( N-1, BIGNUM, X, INCX )
            BETA = BETA*BIGNUM
            ALPHI = ALPHI*BIGNUM
            ALPHR = ALPHR*BIGNUM
            IF( (ABS( BETA ).LT.SMLNUM) .AND. (KNT .LT. 20) )
     $         GO TO 10
*
*           New BETA is at most 1, at least SMLNUM
*
            XNORM = SCNRM2( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI )
            BETA = SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         SAVEALPHA = ALPHA
         ALPHA = ALPHA + BETA
         IF( BETA.LT.ZERO ) THEN
            BETA = -BETA
            TAU = -ALPHA / BETA
         ELSE
            ALPHR = ALPHI * (ALPHI/REAL( ALPHA ))
            ALPHR = ALPHR + XNORM * (XNORM/REAL( ALPHA ))
            TAU = CMPLX( ALPHR/BETA, -ALPHI/BETA )
            ALPHA = CMPLX( -ALPHR, ALPHI )
         END IF
         ALPHA = CLADIV( CMPLX( ONE ), ALPHA )
*
         IF ( ABS(TAU).LE.SMLNUM ) THEN
*
*           In the case where the computed TAU ends up being a denormalized number,
*           it loses relative accuracy. This is a BIG problem. Solution: flush TAU
*           to ZERO (or TWO or whatever makes a nonnegative real number for BETA).
*
*           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
*           (Thanks Pat. Thanks MathWorks.)
*
            ALPHR = REAL( SAVEALPHA )
            ALPHI = AIMAG( SAVEALPHA )
            IF( ALPHI.EQ.ZERO ) THEN
               IF( ALPHR.GE.ZERO ) THEN
                  TAU = ZERO
               ELSE
                  TAU = TWO
                  DO J = 1, N-1
                     X( 1 + (J-1)*INCX ) = ZERO
                  END DO
                  BETA = REAL( -SAVEALPHA )
               END IF
            ELSE
               XNORM = SLAPY2( ALPHR, ALPHI )
               TAU = CMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM )
               DO J = 1, N-1
                  X( 1 + (J-1)*INCX ) = ZERO
               END DO
               BETA = XNORM
            END IF
*
         ELSE
*
*           This is the general case.
*
            CALL CSCAL( N-1, ALPHA, X, INCX )
*
         END IF
*
*        If BETA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SMLNUM
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of CLARFGP
*
      END