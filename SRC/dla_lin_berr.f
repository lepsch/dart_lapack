      SUBROUTINE DLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                N, NZ, NRHS;
      // ..
      // .. Array Arguments ..
      double             AYB( N, NRHS ), BERR( NRHS );
      double             RES( N, NRHS );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      double             TMP;
      int                I, J;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      // EXTERNAL DLAMCH
      double             DLAMCH;
      double             SAFE1;
      // ..
      // .. Executable Statements ..
*
      // Adding SAFE1 to the numerator guards against spuriously zero
      // residuals.  A similar safeguard is in the SLA_yyAMV routine used
     t // o compute AYB.
*
      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (NZ+1)*SAFE1

      DO J = 1, NRHS
         BERR(J) = 0.0D+0
         DO I = 1, N
            IF (AYB(I,J) .NE. 0.0D+0) THEN
               TMP = (SAFE1+ABS(RES(I,J)))/AYB(I,J)
               BERR(J) = MAX( BERR(J), TMP )
            END IF
*
      // If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
     t // he true residual also must be exactly 0.0.
*
         END DO
      END DO
*
      // End of DLA_LIN_BERR
*
      END
