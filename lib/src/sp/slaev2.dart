      void slaev2(A, B, C, RT1, RT2, CS1, SN1 ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               A, B, C, CS1, RT1, RT2, SN1;
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      double               TWO;
      const              TWO = 2.0 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      double               HALF;
      const              HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      int                SGN1, SGN2;
      double               AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      // Compute the eigenvalues

      SM = A + C;
      DF = A - C;
      ADF = ( DF ).abs();
      TB = B + B;
      AB = ( TB ).abs();
      if ( ( A ).abs() > ( C ).abs() ) {
         ACMX = A;
         ACMN = C;
      } else {
         ACMX = C;
         ACMN = A;
      }
      if ( ADF > AB ) {
         RT = ADF*sqrt( ONE+( AB / ADF )**2 );
      } else if ( ADF < AB ) {
         RT = AB*sqrt( ONE+( ADF / AB )**2 );
      } else {

         // Includes case AB=ADF=0

         RT = AB*sqrt( TWO );
      }
      if ( SM < ZERO ) {
         RT1 = HALF*( SM-RT );
         SGN1 = -1;

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
      } else if ( SM > ZERO ) {
         RT1 = HALF*( SM+RT );
         SGN1 = 1;

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
      } else {

         // Includes case RT1 = RT2 = 0

         RT1 = HALF*RT;
         RT2 = -HALF*RT;
         SGN1 = 1;
      }

      // Compute the eigenvector

      if ( DF >= ZERO ) {
         CS = DF + RT;
         SGN2 = 1;
      } else {
         CS = DF - RT;
         SGN2 = -1;
      }
      ACS = ( CS ).abs();
      if ( ACS > AB ) {
         CT = -TB / CS;
         SN1 = ONE / sqrt( ONE+CT*CT );
         CS1 = CT*SN1;
      } else {
         if ( AB == ZERO ) {
            CS1 = ONE;
            SN1 = ZERO;
         } else {
            TN = -CS / TB;
            CS1 = ONE / sqrt( ONE+TN*TN );
            SN1 = TN*CS1;
         }
      }
      if ( SGN1 == SGN2 ) {
         TN = CS1;
         CS1 = -SN1;
         SN1 = TN;
      }
      return;
      }