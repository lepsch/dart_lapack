      REAL             FUNCTION SROUNDUP_LWORK( LWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int              LWORK;
      // ..

* =====================================================================
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC EPSILON, REAL, INT
      // ..
      // .. Executable Statements ..
      // ..
      SROUNDUP_LWORK = REAL( LWORK )

      if ( INT( SROUNDUP_LWORK ) < LWORK ) {
          // Force round up of LWORK
          SROUNDUP_LWORK = SROUNDUP_LWORK * ( 1.0E+0 + EPSILON(0.0E+0) )
      }

      RETURN

      // End of SROUNDUP_LWORK

      }
