      SUBROUTINE DLAE2( A, B, C, RT1, RT2 );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, RT1, RT2;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      double             TWO;
      const              TWO = 2.0 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      double             HALF;
      const              HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      double             AB, ACMN, ACMX, ADF, DF, RT, SM, TB;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      // Compute the eigenvalues

      SM = A + C;
      DF = A - C;
      ADF = ABS( DF );
      TB = B + B;
      AB = ABS( TB );
      if ( ABS( A ) > ABS( C ) ) {
         ACMX = A;
         ACMN = C;
      } else {
         ACMX = C;
         ACMN = A;
      }
      if ( ADF > AB ) {
         RT = ADF*SQRT( ONE+( AB / ADF )**2 );
      } else if ( ADF < AB ) {
         RT = AB*SQRT( ONE+( ADF / AB )**2 );
      } else {

         // Includes case AB=ADF=0

         RT = AB*SQRT( TWO );
      }
      if ( SM < ZERO ) {
         RT1 = HALF*( SM-RT );

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
      } else if ( SM > ZERO ) {
         RT1 = HALF*( SM+RT );

         // Order of execution important.
         // To get fully accurate smaller eigenvalue,
         // next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B;
      } else {

         // Includes case RT1 = RT2 = 0

         RT1 = HALF*RT;
         RT2 = -HALF*RT;
      }
      return;

      // End of DLAE2

      }
