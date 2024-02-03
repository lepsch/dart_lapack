      SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               A, B, C, CS1, RT1, RT2, SN1
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               TWO
      const              TWO = 2.0E0 ;
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               HALF
      const              HALF = 0.5E0 ;
      // ..
      // .. Local Scalars ..
      int                SGN1, SGN2;
      REAL               AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      // Compute the eigenvalues

      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      } else {
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      } else {

         // Includes case AB=ADF=0

         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      } else {

         // Includes case RT1 = RT2 = 0

         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF

      // Compute the eigenvector

      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      } else {
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      } else {
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         } else {
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN

      // End of SLAEV2

      }
