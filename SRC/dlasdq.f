      void dlasdq(UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               ROTATE;
      int                I, ISUB, IUPLO, J, NP1, SQRE1;
      double             CS, R, SMIN, SN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDSQR, DLARTG, DLASR, DSWAP, XERBLA
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      IUPLO = 0;
      if( LSAME( UPLO, 'U' ) ) IUPLO = 1;
      IF( LSAME( UPLO, 'L' ) ) IUPLO = 2;
      if ( IUPLO == 0 ) {
         INFO = -1;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NCVT < 0 ) {
         INFO = -4;
      } else if ( NRU < 0 ) {
         INFO = -5;
      } else if ( NCC < 0 ) {
         INFO = -6;
      } else if ( ( NCVT == 0 && LDVT < 1 ) || ( NCVT > 0 && LDVT < max( 1, N ) ) ) {
         INFO = -10;
      } else if ( LDU < max( 1, NRU ) ) {
         INFO = -12;
      } else if ( ( NCC == 0 && LDC < 1 ) || ( NCC > 0 && LDC < max( 1, N ) ) ) {
         INFO = -14;
      }
      if ( INFO != 0 ) {
         xerbla('DLASDQ', -INFO );
         return;
      }
      if (N == 0) return;

      // ROTATE is true if any singular vectors desired, false otherwise

      ROTATE = ( NCVT > 0 ) || ( NRU > 0 ) || ( NCC > 0 );
      NP1 = N + 1;
      SQRE1 = SQRE;

      // If matrix non-square upper bidiagonal, rotate to be lower
      // bidiagonal.  The rotations are on the right.

      if ( ( IUPLO == 1 ) && ( SQRE1 == 1 ) ) {
         for (I = 1; I <= N - 1; I++) { // 10
            dlartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R;
            E( I ) = SN*D( I+1 );
            D( I+1 ) = CS*D( I+1 );
            if ( ROTATE ) {
               WORK( I ) = CS;
               WORK( N+I ) = SN;
            }
         } // 10
         dlartg(D( N ), E( N ), CS, SN, R );
         D( N ) = R;
         E( N ) = ZERO;
         if ( ROTATE ) {
            WORK( N ) = CS;
            WORK( N+N ) = SN;
         }
         IUPLO = 2;
         SQRE1 = 0;

         // Update singular vectors if desired.

         if (NCVT > 0) dlasr( 'L', 'V', 'F', NP1, NCVT, WORK( 1 ), WORK( NP1 ), VT, LDVT );
      }

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left.

      if ( IUPLO == 2 ) {
         for (I = 1; I <= N - 1; I++) { // 20
            dlartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R;
            E( I ) = SN*D( I+1 );
            D( I+1 ) = CS*D( I+1 );
            if ( ROTATE ) {
               WORK( I ) = CS;
               WORK( N+I ) = SN;
            }
         } // 20

         // If matrix (N+1)-by-N lower bidiagonal, one additional
         // rotation is needed.

         if ( SQRE1 == 1 ) {
            dlartg(D( N ), E( N ), CS, SN, R );
            D( N ) = R;
            if ( ROTATE ) {
               WORK( N ) = CS;
               WORK( N+N ) = SN;
            }
         }

         // Update singular vectors if desired.

         if ( NRU > 0 ) {
            if ( SQRE1 == 0 ) {
               dlasr('R', 'V', 'F', NRU, N, WORK( 1 ), WORK( NP1 ), U, LDU );
            } else {
               dlasr('R', 'V', 'F', NRU, NP1, WORK( 1 ), WORK( NP1 ), U, LDU );
            }
         }
         if ( NCC > 0 ) {
            if ( SQRE1 == 0 ) {
               dlasr('L', 'V', 'F', N, NCC, WORK( 1 ), WORK( NP1 ), C, LDC );
            } else {
               dlasr('L', 'V', 'F', NP1, NCC, WORK( 1 ), WORK( NP1 ), C, LDC );
            }
         }
      }

      // Call DBDSQR to compute the SVD of the reduced real
      // N-by-N upper bidiagonal matrix.

      dbdsqr('U', N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO );

      // Sort the singular values into ascending order (insertion sort on
      // singular values, but only one transposition per singular vector)

      for (I = 1; I <= N; I++) { // 40

         // Scan for smallest D(I).

         ISUB = I;
         SMIN = D( I );
         for (J = I + 1; J <= N; J++) { // 30
            if ( D( J ) < SMIN ) {
               ISUB = J;
               SMIN = D( J );
            }
         } // 30
         if ( ISUB != I ) {

            // Swap singular values and vectors.

            D( ISUB ) = D( I );
            D( I ) = SMIN;
            if (NCVT > 0) dswap( NCVT, VT( ISUB, 1 ), LDVT, VT( I, 1 ), LDVT );
            if( NRU > 0 ) dswap( NRU, U( 1, ISUB ), 1, U( 1, I ), 1 );
            IF( NCC > 0 ) dswap( NCC, C( ISUB, 1 ), LDC, C( I, 1 ), LDC );
         }
      } // 40

      return;
      }
