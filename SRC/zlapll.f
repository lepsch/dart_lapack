      void zlapll(N, X, INCX, Y, INCY, SSMIN ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      double             SSMIN;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      COMPLEX*16         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      double             SSMAX;
      COMPLEX*16         A11, A12, A22, C, TAU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG
      // ..
      // .. External Functions ..
      COMPLEX*16         ZDOTC;
      // EXTERNAL ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAS2, ZAXPY, ZLARFG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 1 ) {
         SSMIN = ZERO;
         return;
      }

      // Compute the QR factorization of the N-by-2 matrix ( X Y )

      zlarfg(N, X( 1 ), X( 1+INCX ), INCX, TAU );
      A11 = X( 1 );
      X( 1 ) = CONE;

      C = -DCONJG( TAU )*ZDOTC( N, X, INCX, Y, INCY );
      zaxpy(N, C, X, INCX, Y, INCY );

      zlarfg(N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU );

      A12 = Y( 1 );
      A22 = Y( 1+INCY );

      // Compute the SVD of 2-by-2 Upper triangular matrix.

      dlas2(ABS( A11 ), ABS( A12 ), ABS( A22 ), SSMIN, SSMAX );

      return;

      // End of ZLAPLL

      }
