      SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      COMPLEX            ALPHA, TAU
      // ..
      // .. Array Arguments ..
      COMPLEX            X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                J, KNT;
      REAL               ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
      // ..
      // .. External Functions ..
      REAL               SCNRM2, SLAMCH, SLAPY3
      COMPLEX            CLADIV
      // EXTERNAL SCNRM2, SLAMCH, SLAPY3, CLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, REAL, SIGN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSSCAL
      // ..
      // .. Executable Statements ..

      if ( N <= 0 ) {
         TAU = ZERO
         RETURN
      }

      XNORM = SCNRM2( N-1, X, INCX )
      ALPHR = REAL( ALPHA )
      ALPHI = AIMAG( ALPHA )

      if ( XNORM == ZERO && ALPHI == ZERO ) {

         // H  =  I

         TAU = ZERO
      } else {

         // general case

         BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN

         KNT = 0
         if ( ABS( BETA ) < SAFMIN ) {

            // XNORM, BETA may be inaccurate; scale X and recompute them

            } // 10
            KNT = KNT + 1
            csscal(N-1, RSAFMN, X, INCX );
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( (ABS( BETA ) < SAFMIN) && (KNT < 20) ) GO TO 10

            // New BETA is at most 1, at least SAFMIN

            XNORM = SCNRM2( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI )
            BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         }
         TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
         cscal(N-1, ALPHA, X, INCX );

         // If ALPHA is subnormal, it may lose relative accuracy

         for (J = 1; J <= KNT; J++) { // 20
            BETA = BETA*SAFMIN
         } // 20
         ALPHA = BETA
      }

      RETURN

      // End of CLARFG

      }
