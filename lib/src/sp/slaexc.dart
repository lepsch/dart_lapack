      void slaexc(WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ;
      int                INFO, J1, LDQ, LDT, N, N1, N2;
      // ..
      // .. Array Arguments ..
      double               Q( LDQ, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               TEN;
      const              TEN = 1.0e+1 ;
      int                LDD, LDX;
      const              LDD = 4, LDX = 2 ;
      // ..
      // .. Local Scalars ..
      int                IERR, J2, J3, J4, K, ND;
      double               CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22, T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2, WR1, WR2, XNORM;
      // ..
      // .. Local Arrays ..
      double               D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ), X( LDX, 2 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLANV2, SLARFG, SLARFX, SLARTG, SLASY2, SROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // Quick return if possible

      if (N == 0 || N1 == 0 || N2 == 0) return;
      IF( J1+N1 > N ) return;

      J2 = J1 + 1;
      J3 = J1 + 2;
      J4 = J1 + 3;

      if ( N1 == 1 && N2 == 1 ) {

         // Swap two 1-by-1 blocks.

         T11 = T( J1, J1 );
         T22 = T( J2, J2 );

         // Determine the transformation to perform the interchange.

         slartg(T( J1, J2 ), T22-T11, CS, SN, TEMP );

         // Apply transformation to the matrix T.

         if (J3 <= N) srot( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS, SN );
         srot(J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN );

         T[J1, J1] = T22;
         T[J2, J2] = T11;

         if ( WANTQ ) {

            // Accumulate transformation in the matrix Q.

            srot(N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN );
         }

      } else {

         // Swapping involves at least one 2-by-2 block.

         // Copy the diagonal block of order N1+N2 to the local array D
         // and compute its norm.

         ND = N1 + N2;
         slacpy('Full', ND, ND, T( J1, J1 ), LDT, D, LDD );
         DNORM = SLANGE( 'Max', ND, ND, D, LDD, WORK );

         // Compute machine-dependent threshold for test for accepting
         // swap.

         EPS = SLAMCH( 'P' );
         SMLNUM = SLAMCH( 'S' ) / EPS;
         THRESH = max( TEN*EPS*DNORM, SMLNUM );

         // Solve T11*X - X*T22 = scale*T12 for X.

         slasy2( false , false , -1, N1, N2, D, LDD, D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X, LDX, XNORM, IERR );

         // Swap the adjacent diagonal blocks.

         K = N1 + N1 + N2 - 3;
         GO TO ( 10, 20, 30 )K;

         } // 10

         // N1 = 1, N2 = 2: generate elementary reflector H so that:

         // ( scale, X11, X12 ) H = ( 0, 0, * )

         U[1] = SCALE;
         U[2] = X( 1, 1 );
         U[3] = X( 1, 2 );
         slarfg(3, U( 3 ), U, 1, TAU );
         U[3] = ONE;
         T11 = T( J1, J1 );

         // Perform swap provisionally on diagonal block in D.

         slarfx('L', 3, 3, U, TAU, D, LDD, WORK );
         slarfx('R', 3, 3, U, TAU, D, LDD, WORK );

         // Test whether to reject swap.

         if( max( ( D( 3, 1 ) ).abs(), ( D( 3, 2 ) ).abs(), ( D( 3, 3 )-T11 ) ).abs() > THRESH )GO TO 50;

         // Accept swap: apply transformation to the entire matrix T.

         slarfx('L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK );
         slarfx('R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK );

         T[J3, J1] = ZERO;
         T[J3, J2] = ZERO;
         T[J3, J3] = T11;

         if ( WANTQ ) {

            // Accumulate transformation in the matrix Q.

            slarfx('R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK );
         }
         GO TO 40;

         } // 20

         // N1 = 2, N2 = 1: generate elementary reflector H so that:

         // H (  -X11 ) = ( * )
           // (  -X21 ) = ( 0 )
           // ( scale ) = ( 0 )

         U[1] = -X( 1, 1 );
         U[2] = -X( 2, 1 );
         U[3] = SCALE;
         slarfg(3, U( 1 ), U( 2 ), 1, TAU );
         U[1] = ONE;
         T33 = T( J3, J3 );

         // Perform swap provisionally on diagonal block in D.

         slarfx('L', 3, 3, U, TAU, D, LDD, WORK );
         slarfx('R', 3, 3, U, TAU, D, LDD, WORK );

         // Test whether to reject swap.

         if( max( ( D( 2, 1 ) ).abs(), ( D( 3, 1 ) ).abs(), ( D( 1, 1 )-T33 ) ).abs() > THRESH )GO TO 50;

         // Accept swap: apply transformation to the entire matrix T.

         slarfx('R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK );
         slarfx('L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK );

         T[J1, J1] = T33;
         T[J2, J1] = ZERO;
         T[J3, J1] = ZERO;

         if ( WANTQ ) {

            // Accumulate transformation in the matrix Q.

            slarfx('R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK );
         }
         GO TO 40;

         } // 30

         // N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
         // that:

         // H(2) H(1) (  -X11  -X12 ) = (  *  * )
                   // (  -X21  -X22 )   (  0  * )
                   // ( scale    0  )   (  0  0 )
                   // (    0  scale )   (  0  0 )

         U1[1] = -X( 1, 1 );
         U1[2] = -X( 2, 1 );
         U1[3] = SCALE;
         slarfg(3, U1( 1 ), U1( 2 ), 1, TAU1 );
         U1[1] = ONE;

         TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) );
         U2[1] = -TEMP*U1( 2 ) - X( 2, 2 );
         U2[2] = -TEMP*U1( 3 );
         U2[3] = SCALE;
         slarfg(3, U2( 1 ), U2( 2 ), 1, TAU2 );
         U2[1] = ONE;

         // Perform swap provisionally on diagonal block in D.

         slarfx('L', 3, 4, U1, TAU1, D, LDD, WORK );
         slarfx('R', 4, 3, U1, TAU1, D, LDD, WORK );
         slarfx('L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK );
         slarfx('R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK );

         // Test whether to reject swap.

         if( max( ( D( 3, 1 ) ).abs(), ( D( 3, 2 ) ).abs(), ( D( 4, 1 ) ).abs(), ( D( 4, 2 ) ) ).abs() > THRESH )GO TO 50;

         // Accept swap: apply transformation to the entire matrix T.

         slarfx('L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK );
         slarfx('R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK );
         slarfx('L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK );
         slarfx('R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK );

         T[J3, J1] = ZERO;
         T[J3, J2] = ZERO;
         T[J4, J1] = ZERO;
         T[J4, J2] = ZERO;

         if ( WANTQ ) {

            // Accumulate transformation in the matrix Q.

            slarfx('R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK );
            slarfx('R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK );
         }

         } // 40

         if ( N2 == 2 ) {

            // Standardize new 2-by-2 block T11

            slanv2(T( J1, J1 ), T( J1, J2 ), T( J2, J1 ), T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN );
            srot(N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT, CS, SN );
            srot(J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN );
            if (WANTQ) srot( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN );
         }

         if ( N1 == 2 ) {

            // Standardize new 2-by-2 block T22

            J3 = J1 + N2;
            J4 = J3 + 1;
            slanv2(T( J3, J3 ), T( J3, J4 ), T( J4, J3 ), T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )             IF( J3+2 <= N ) CALL SROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ), LDT, CS, SN );
            srot(J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN );
            if (WANTQ) srot( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN );
         }

      }
      return;

      // Exit with INFO = 1 if swap was rejected.

   50 INFO = 1;
      return;
      }
