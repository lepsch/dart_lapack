      SUBROUTINE ZLAPLL( N, X, INCX, Y, INCY, SSMIN )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      double             SSMIN;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      double             SSMAX;
      COMPLEX*16         A11, A12, A22, C, TAU
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG
      // ..
      // .. External Functions ..
      COMPLEX*16         ZDOTC
      // EXTERNAL ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAS2, ZAXPY, ZLARFG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.1 ) THEN
         SSMIN = ZERO
         RETURN
      END IF

      // Compute the QR factorization of the N-by-2 matrix ( X Y )

      CALL ZLARFG( N, X( 1 ), X( 1+INCX ), INCX, TAU )
      A11 = X( 1 )
      X( 1 ) = CONE

      C = -DCONJG( TAU )*ZDOTC( N, X, INCX, Y, INCY )
      CALL ZAXPY( N, C, X, INCX, Y, INCY )

      CALL ZLARFG( N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU )

      A12 = Y( 1 )
      A22 = Y( 1+INCY )

      // Compute the SVD of 2-by-2 Upper triangular matrix.

      CALL DLAS2( ABS( A11 ), ABS( A12 ), ABS( A22 ), SSMIN, SSMAX )

      RETURN

      // End of ZLAPLL

      }
