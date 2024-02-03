      SUBROUTINE ZLATSY( UPLO, N, X, LDX, ISEED );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDX, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( * );
      COMPLEX*16         X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         EYE;
      const              EYE = ( 0.0, 1.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, N5;
      double             ALPHA, ALPHA3, BETA;
      COMPLEX*16         A, B, C, R;
      // ..
      // .. External Functions ..
      COMPLEX*16         ZLARND;
      // EXTERNAL ZLARND
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      ALPHA = ( 1.0+SQRT( 17.0 ) ) / 8.0;
      BETA = ALPHA - 1.0 / 1000.0;
      ALPHA3 = ALPHA*ALPHA*ALPHA;

      // UPLO = 'U':  Upper triangular storage

      if ( UPLO == 'U' ) {

         // Fill the upper triangle of the matrix with zeros.

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               X( I, J ) = 0.0;
            } // 10
         } // 20
         N5 = N / 5;
         N5 = N - 5*N5 + 1;

         DO 30 I = N, N5, -5;
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( I, I ) = A;
            X( I-2, I ) = B;
            X( I-2, I-1 ) = R;
            X( I-2, I-2 ) = C;
            X( I-1, I-1 ) = ZLARND( 2, ISEED );
            X( I-3, I-3 ) = ZLARND( 2, ISEED );
            X( I-4, I-4 ) = ZLARND( 2, ISEED );
            if ( ABS( X( I-3, I-3 ) ) > ABS( X( I-4, I-4 ) ) ) {
               X( I-4, I-3 ) = 2.0*X( I-3, I-3 );
            } else {
               X( I-4, I-3 ) = 2.0*X( I-4, I-4 );
            }
         } // 30

         // Clean-up for N not a multiple of 5.

         I = N5 - 1;
         if ( I > 2 ) {
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( I, I ) = A;
            X( I-2, I ) = B;
            X( I-2, I-1 ) = R;
            X( I-2, I-2 ) = C;
            X( I-1, I-1 ) = ZLARND( 2, ISEED );
            I = I - 3;
         }
         if ( I > 1 ) {
            X( I, I ) = ZLARND( 2, ISEED );
            X( I-1, I-1 ) = ZLARND( 2, ISEED );
            if ( ABS( X( I, I ) ) > ABS( X( I-1, I-1 ) ) ) {
               X( I-1, I ) = 2.0*X( I, I );
            } else {
               X( I-1, I ) = 2.0*X( I-1, I-1 );
            }
            I = I - 2;
         } else if ( I == 1 ) {
            X( I, I ) = ZLARND( 2, ISEED );
            I = I - 1;
         }

      // UPLO = 'L':  Lower triangular storage

      } else {

         // Fill the lower triangle of the matrix with zeros.

         for (J = 1; J <= N; J++) { // 50
            for (I = J; I <= N; I++) { // 40
               X( I, J ) = 0.0;
            } // 40
         } // 50
         N5 = N / 5;
         N5 = N5*5;

         DO 60 I = 1, N5, 5;
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( I, I ) = A;
            X( I+2, I ) = B;
            X( I+2, I+1 ) = R;
            X( I+2, I+2 ) = C;
            X( I+1, I+1 ) = ZLARND( 2, ISEED );
            X( I+3, I+3 ) = ZLARND( 2, ISEED );
            X( I+4, I+4 ) = ZLARND( 2, ISEED );
            if ( ABS( X( I+3, I+3 ) ) > ABS( X( I+4, I+4 ) ) ) {
               X( I+4, I+3 ) = 2.0*X( I+3, I+3 );
            } else {
               X( I+4, I+3 ) = 2.0*X( I+4, I+4 );
            }
         } // 60

         // Clean-up for N not a multiple of 5.

         I = N5 + 1;
         if ( I < N-1 ) {
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( I, I ) = A;
            X( I+2, I ) = B;
            X( I+2, I+1 ) = R;
            X( I+2, I+2 ) = C;
            X( I+1, I+1 ) = ZLARND( 2, ISEED );
            I = I + 3;
         }
         if ( I < N ) {
            X( I, I ) = ZLARND( 2, ISEED );
            X( I+1, I+1 ) = ZLARND( 2, ISEED );
            if ( ABS( X( I, I ) ) > ABS( X( I+1, I+1 ) ) ) {
               X( I+1, I ) = 2.0*X( I, I );
            } else {
               X( I+1, I ) = 2.0*X( I+1, I+1 );
            }
            I = I + 2;
         } else if ( I == N ) {
            X( I, I ) = ZLARND( 2, ISEED );
            I = I + 1;
         }
      }

      RETURN;

      // End of ZLATSY

      }
