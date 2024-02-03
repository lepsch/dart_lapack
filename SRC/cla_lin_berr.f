      SUBROUTINE CLA_LIN_BERR( N, NZ, NRHS, RES, AYB, BERR )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NZ, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               AYB( N, NRHS ), BERR( NRHS )
      COMPLEX            RES( N, NRHS )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      REAL               TMP
      int                I, J;
      COMPLEX            CDUM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, AIMAG, MAX
      // ..
      // .. External Functions ..
      // EXTERNAL SLAMCH
      REAL               SLAMCH
      REAL               SAFE1
      // ..
      // .. Statement Functions ..
      COMPLEX            CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Adding SAFE1 to the numerator guards against spuriously zero
      // residuals.  A similar safeguard is in the CLA_yyAMV routine used
      // to compute AYB.

      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (NZ+1)*SAFE1

      DO J = 1, NRHS
         BERR(J) = 0.0
         DO I = 1, N
            if (AYB(I,J) .NE. 0.0) {
               TMP = (SAFE1 + CABS1(RES(I,J)))/AYB(I,J)
               BERR(J) = MAX( BERR(J), TMP )
            }

      // If AYB is exactly 0.0 (and if computed by CLA_yyAMV), then we know
      // the true residual also must be exactly 0.0.

         END DO
      END DO

      // End of CLA_LIN_BERR

      }
