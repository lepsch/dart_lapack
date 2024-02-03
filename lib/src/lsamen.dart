      bool lsamen(N, CA, CB ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       CA, CB;
      int                N;
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN
      // ..
      // .. Executable Statements ..

      LSAMEN = false;
      if( LEN( CA ) < N || LEN( CB ) < N ) GO TO 20;

      // Do for each character in the two strings.

      for (I = 1; I <= N; I++) { // 10

         // Test if the characters are equal using LSAME.

         if( !LSAME( CA( I: I ), CB( I: I ) ) ) GO TO 20;

      } // 10
      LSAMEN = true;

      } // 20
      return;
      }
