      SUBROUTINE DLAPLL( N, X, INCX, Y, INCY, SSMIN )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      double             SSMIN;
      // ..
      // .. Array Arguments ..
      double             X( * ), Y( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             A11, A12, A22, C, SSMAX, TAU;
      // ..
      // .. External Functions ..
      double             DDOT;
      // EXTERNAL DDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLARFG, DLAS2
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.1 ) {
         SSMIN = ZERO
         RETURN
      }

      // Compute the QR factorization of the N-by-2 matrix ( X Y )

      dlarfg(N, X( 1 ), X( 1+INCX ), INCX, TAU );
      A11 = X( 1 )
      X( 1 ) = ONE

      C = -TAU*DDOT( N, X, INCX, Y, INCY )
      daxpy(N, C, X, INCX, Y, INCY );

      dlarfg(N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU );

      A12 = Y( 1 )
      A22 = Y( 1+INCY )

      // Compute the SVD of 2-by-2 Upper triangular matrix.

      dlas2(A11, A12, A22, SSMIN, SSMAX );

      RETURN

      // End of DLAPLL

      }
