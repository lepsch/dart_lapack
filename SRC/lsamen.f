      bool             FUNCTION LSAMEN( N, CA, CB );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       CA, CB;
      int                N;
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN
      // ..
      // .. Executable Statements ..

      LSAMEN = false;
      IF( LEN( CA ).LT.N || LEN( CB ).LT.N ) GO TO 20

      // Do for each character in the two strings.

      for (I = 1; I <= N; I++) { // 10

         // Test if the characters are equal using LSAME.

         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) ) GO TO 20

      } // 10
      LSAMEN = true;

      } // 20
      RETURN

      // End of LSAMEN

      }
