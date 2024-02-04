      void clapll(N, X, INCX, Y, INCY, SSMIN ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      REAL               SSMIN;
      // ..
      // .. Array Arguments ..
      COMPLEX            X( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      REAL               SSMAX;
      COMPLEX            A11, A12, A22, C, TAU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG
      // ..
      // .. External Functions ..
      //- COMPLEX            CDOTC;
      // EXTERNAL CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CLARFG, SLAS2
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 1 ) {
         SSMIN = ZERO;
         return;
      }

      // Compute the QR factorization of the N-by-2 matrix ( X Y )

      clarfg(N, X( 1 ), X( 1+INCX ), INCX, TAU );
      A11 = X( 1 );
      X[1] = CONE;

      C = -CONJG( TAU )*CDOTC( N, X, INCX, Y, INCY );
      caxpy(N, C, X, INCX, Y, INCY );

      clarfg(N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU );

      A12 = Y( 1 );
      A22 = Y( 1+INCY );

      // Compute the SVD of 2-by-2 Upper triangular matrix.

      slas2(( A11 ).abs(), ( A12 ).abs(), ( A22 ).abs(), SSMIN, SSMAX );

      return;
      }
