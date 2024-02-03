      SUBROUTINE ZLATSP( UPLO, N, X, ISEED );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      int                ISEED( * );
      COMPLEX*16         X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         EYE;
      const              EYE = ( 0.0, 1.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, JJ, N5;
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

      // Fill the matrix with zeros.

      for (J = 1; J <= N*( N+1 ) / 2; J++) { // 10
         X( J ) = 0.0;
      } // 10

      // UPLO = 'U':  Upper triangular storage

      if ( UPLO == 'U' ) {
         N5 = N / 5;
         N5 = N - 5*N5 + 1;

         JJ = N*( N+1 ) / 2;
         DO 20 J = N, N5, -5;
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( JJ ) = A;
            X( JJ-2 ) = B;
            JJ = JJ - J;
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ-1 ) = R;
            JJ = JJ - ( J-1 );
            X( JJ ) = C;
            JJ = JJ - ( J-2 );
            X( JJ ) = ZLARND( 2, ISEED );
            JJ = JJ - ( J-3 );
            X( JJ ) = ZLARND( 2, ISEED );
            if ( ABS( X( JJ+( J-3 ) ) ) > ABS( X( JJ ) ) ) {
               X( JJ+( J-4 ) ) = 2.0*X( JJ+( J-3 ) );
            } else {
               X( JJ+( J-4 ) ) = 2.0*X( JJ );
            }
            JJ = JJ - ( J-4 );
         } // 20

         // Clean-up for N not a multiple of 5.

         J = N5 - 1;
         if ( J > 2 ) {
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( JJ ) = A;
            X( JJ-2 ) = B;
            JJ = JJ - J;
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ-1 ) = R;
            JJ = JJ - ( J-1 );
            X( JJ ) = C;
            JJ = JJ - ( J-2 );
            J = J - 3;
         }
         if ( J > 1 ) {
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ-J ) = ZLARND( 2, ISEED );
            if ( ABS( X( JJ ) ) > ABS( X( JJ-J ) ) ) {
               X( JJ-1 ) = 2.0*X( JJ );
            } else {
               X( JJ-1 ) = 2.0*X( JJ-J );
            }
            JJ = JJ - J - ( J-1 );
            J = J - 2;
         } else if ( J == 1 ) {
            X( JJ ) = ZLARND( 2, ISEED );
            J = J - 1;
         }

      // UPLO = 'L':  Lower triangular storage

      } else {
         N5 = N / 5;
         N5 = N5*5;

         JJ = 1;
         DO 30 J = 1, N5, 5;
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( JJ ) = A;
            X( JJ+2 ) = B;
            JJ = JJ + ( N-J+1 );
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ+1 ) = R;
            JJ = JJ + ( N-J );
            X( JJ ) = C;
            JJ = JJ + ( N-J-1 );
            X( JJ ) = ZLARND( 2, ISEED );
            JJ = JJ + ( N-J-2 );
            X( JJ ) = ZLARND( 2, ISEED );
            if ( ABS( X( JJ-( N-J-2 ) ) ) > ABS( X( JJ ) ) ) {
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ-( N-J-2 ) );
            } else {
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ );
            }
            JJ = JJ + ( N-J-3 );
         } // 30

         // Clean-up for N not a multiple of 5.

         J = N5 + 1;
         if ( J < N-1 ) {
            A = ALPHA3*ZLARND( 5, ISEED );
            B = ZLARND( 5, ISEED ) / ALPHA;
            C = A - 2.0*B*EYE;
            R = C / BETA;
            X( JJ ) = A;
            X( JJ+2 ) = B;
            JJ = JJ + ( N-J+1 );
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ+1 ) = R;
            JJ = JJ + ( N-J );
            X( JJ ) = C;
            JJ = JJ + ( N-J-1 );
            J = J + 3;
         }
         if ( J < N ) {
            X( JJ ) = ZLARND( 2, ISEED );
            X( JJ+( N-J+1 ) ) = ZLARND( 2, ISEED );
            if ( ABS( X( JJ ) ) > ABS( X( JJ+( N-J+1 ) ) ) ) {
               X( JJ+1 ) = 2.0*X( JJ );
            } else {
               X( JJ+1 ) = 2.0*X( JJ+( N-J+1 ) );
            }
            JJ = JJ + ( N-J+1 ) + ( N-J );
            J = J + 2;
         } else if ( J == N ) {
            X( JJ ) = ZLARND( 2, ISEED );
            JJ = JJ + ( N-J+1 );
            J = J + 1;
         }
      }

      RETURN;

      // End of ZLATSP

      }
