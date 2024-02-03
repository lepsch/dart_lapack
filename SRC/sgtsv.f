      SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('SGTSV ', -INFO );
         RETURN;
      }

      if (N == 0) RETURN;

      if ( NRHS == 1 ) {
         for (I = 1; I <= N - 2; I++) { // 10
            if ( ABS( D( I ) ) >= ABS( DL( I ) ) ) {

               // No row interchange required

               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D( I+1 ) = D( I+1 ) - FACT*DU( I );
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 );
               } else {
                  INFO = I;
                  RETURN;
               }
               DL( I ) = ZERO;
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I );
               D( I ) = DL( I );
               TEMP = D( I+1 );
               D( I+1 ) = DU( I ) - FACT*TEMP;
               DL( I ) = DU( I+1 );
               DU( I+1 ) = -FACT*DL( I );
               DU( I ) = TEMP;
               TEMP = B( I, 1 );
               B( I, 1 ) = B( I+1, 1 );
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 );
            }
         } // 10
         if ( N > 1 ) {
            I = N - 1;
            if ( ABS( D( I ) ) >= ABS( DL( I ) ) ) {
               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D( I+1 ) = D( I+1 ) - FACT*DU( I );
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 );
               } else {
                  INFO = I;
                  RETURN;
               }
            } else {
               FACT = D( I ) / DL( I );
               D( I ) = DL( I );
               TEMP = D( I+1 );
               D( I+1 ) = DU( I ) - FACT*TEMP;
               DU( I ) = TEMP;
               TEMP = B( I, 1 );
               B( I, 1 ) = B( I+1, 1 );
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 );
            }
         }
         if ( D( N ) == ZERO ) {
            INFO = N;
            RETURN;
         }
      } else {
         for (I = 1; I <= N - 2; I++) { // 40
            if ( ABS( D( I ) ) >= ABS( DL( I ) ) ) {

               // No row interchange required

               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D( I+1 ) = D( I+1 ) - FACT*DU( I );
                  for (J = 1; J <= NRHS; J++) { // 20
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J );
                  } // 20
               } else {
                  INFO = I;
                  RETURN;
               }
               DL( I ) = ZERO;
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I );
               D( I ) = DL( I );
               TEMP = D( I+1 );
               D( I+1 ) = DU( I ) - FACT*TEMP;
               DL( I ) = DU( I+1 );
               DU( I+1 ) = -FACT*DL( I );
               DU( I ) = TEMP;
               for (J = 1; J <= NRHS; J++) { // 30
                  TEMP = B( I, J );
                  B( I, J ) = B( I+1, J );
                  B( I+1, J ) = TEMP - FACT*B( I+1, J );
               } // 30
            }
         } // 40
         if ( N > 1 ) {
            I = N - 1;
            if ( ABS( D( I ) ) >= ABS( DL( I ) ) ) {
               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D( I+1 ) = D( I+1 ) - FACT*DU( I );
                  for (J = 1; J <= NRHS; J++) { // 50
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J );
                  } // 50
               } else {
                  INFO = I;
                  RETURN;
               }
            } else {
               FACT = D( I ) / DL( I );
               D( I ) = DL( I );
               TEMP = D( I+1 );
               D( I+1 ) = DU( I ) - FACT*TEMP;
               DU( I ) = TEMP;
               for (J = 1; J <= NRHS; J++) { // 60
                  TEMP = B( I, J );
                  B( I, J ) = B( I+1, J );
                  B( I+1, J ) = TEMP - FACT*B( I+1, J );
               } // 60
            }
         }
         if ( D( N ) == ZERO ) {
            INFO = N;
            RETURN;
         }
      }

      // Back solve with the matrix U from the factorization.

      if ( NRHS <= 2 ) {
         J = 1;
         } // 70
         B( N, J ) = B( N, J ) / D( N );
         if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
         DO 80 I = N - 2, 1, -1;
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I );
         } // 80
         if ( J < NRHS ) {
            J = J + 1;
            GO TO 70;
         }
      } else {
         for (J = 1; J <= NRHS; J++) { // 100
            B( N, J ) = B( N, J ) / D( N );
            if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
            DO 90 I = N - 2, 1, -1;
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I );
            } // 90
         } // 100
      }

      RETURN;

      // End of SGTSV

      }
