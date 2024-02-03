      SUBROUTINE DTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, N1, N2, WORK, LWORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================
*  Replaced various illegal calls to DCOPY by calls to DLASET, or by DO
*  loops. Sven Hammarling, 1/5/02.

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             TWENTY;
      const              TWENTY = 2.0e+01 ;
      int                LDST;
      const              LDST = 4 ;
      bool               WANDS;
      const              WANDS = true ;
      // ..
      // .. Local Scalars ..
      bool               STRONG, WEAK;
      int                I, IDUM, LINFO, M;
      double             BQRA21, BRQA21, DDUM, DNORMA, DNORMB, DSCALE, DSUM, EPS, F, G, SA, SB, SCALE, SMLNUM, THRESHA, THRESHB;
      // ..
      // .. Local Arrays ..
      int                IWORK( LDST + 2 );
      double             AI( 2 ), AR( 2 ), BE( 2 ), IR( LDST, LDST ), IRCOP( LDST, LDST ), LI( LDST, LDST ), LICOP( LDST, LDST ), S( LDST, LDST ), SCPY( LDST, LDST ), T( LDST, LDST ), TAUL( LDST ), TAUR( LDST ), TCPY( LDST, LDST );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DGEQR2, DGERQ2, DLACPY, DLAGV2, DLARTG, DLASET, DLASSQ, DORG2R, DORGR2, DORM2R, DORMR2, DROT, DSCAL, DTGSY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if (N <= 1 || N1 <= 0 || N2 <= 0) RETURN       IF( N1 > N || ( J1+N1 ) > N ) RETURN;
      M = N1 + N2
      if ( LWORK < MAX( 1, N*M, M*M*2 ) ) {
         INFO = -16
         WORK( 1 ) = MAX( 1, N*M, M*M*2 )
         RETURN
      }

      WEAK = false;
      STRONG = false;

      // Make a local copy of selected block

      dlaset('Full', LDST, LDST, ZERO, ZERO, LI, LDST );
      dlaset('Full', LDST, LDST, ZERO, ZERO, IR, LDST );
      dlacpy('Full', M, M, A( J1, J1 ), LDA, S, LDST );
      dlacpy('Full', M, M, B( J1, J1 ), LDB, T, LDST );

      // Compute threshold for testing acceptance of swapping.

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      DSCALE = ZERO
      DSUM = ONE
      dlacpy('Full', M, M, S, LDST, WORK, M );
      dlassq(M*M, WORK, 1, DSCALE, DSUM );
      DNORMA = DSCALE*SQRT( DSUM )
      DSCALE = ZERO
      DSUM = ONE
      dlacpy('Full', M, M, T, LDST, WORK, M );
      dlassq(M*M, WORK, 1, DSCALE, DSUM );
      DNORMB = DSCALE*SQRT( DSUM )

      // THRES has been changed from
         // THRESH = MAX( TEN*EPS*SA, SMLNUM )
      // to
         // THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
      // on 04/01/10.
      // "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
      // Jim Demmel and Guillaume Revy. See forum post 1783.

      THRESHA = MAX( TWENTY*EPS*DNORMA, SMLNUM )
      THRESHB = MAX( TWENTY*EPS*DNORMB, SMLNUM )

      if ( M == 2 ) {

         // CASE 1: Swap 1-by-1 and 1-by-1 blocks.

         // Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
         // using Givens rotations and perform the swap tentatively.

         F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
         G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
         SA = ABS( S( 2, 2 ) ) * ABS( T( 1, 1 ) )
         SB = ABS( S( 1, 1 ) ) * ABS( T( 2, 2 ) )
         dlartg(F, G, IR( 1, 2 ), IR( 1, 1 ), DDUM );
         IR( 2, 1 ) = -IR( 1, 2 )
         IR( 2, 2 ) = IR( 1, 1 )
         drot(2, S( 1, 1 ), 1, S( 1, 2 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         drot(2, T( 1, 1 ), 1, T( 1, 2 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         if ( SA >= SB ) {
            dlartg(S( 1, 1 ), S( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), DDUM );
         } else {
            dlartg(T( 1, 1 ), T( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), DDUM );
         }
         drot(2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, LI( 1, 1 ), LI( 2, 1 ) );
         drot(2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, LI( 1, 1 ), LI( 2, 1 ) );
         LI( 2, 2 ) = LI( 1, 1 )
         LI( 1, 2 ) = -LI( 2, 1 )

         // Weak stability test: |S21| <= O(EPS F-norm((A)))
                            // and  |T21| <= O(EPS F-norm((B)))

         WEAK = ABS( S( 2, 1 ) ) <= THRESHA && ABS( T( 2, 1 ) ) <= THRESHB          IF( .NOT.WEAK ) GO TO 70

         if ( WANDS ) {

            // Strong stability test:
                // F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
                // and
                // F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))

            dlacpy('Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), M );
            dgemm('N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
            dgemm('N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            dlassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SA = DSCALE*SQRT( DSUM )

            dlacpy('Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), M );
            dgemm('N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
            dgemm('N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            dlassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SB = DSCALE*SQRT( DSUM )
            STRONG = SA <= THRESHA && SB <= THRESHB
            if (.NOT.STRONG) GO TO 70;
         }

         // Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
                // (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).

         drot(J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         drot(J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         drot(N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, LI( 1, 1 ), LI( 2, 1 ) );
         drot(N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, LI( 1, 1 ), LI( 2, 1 ) );

         // Set  N1-by-N2 (2,1) - blocks to ZERO.

         A( J1+1, J1 ) = ZERO
         B( J1+1, J1 ) = ZERO

         // Accumulate transformations into Q and Z if requested.

         if (WANTZ) CALL DROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) )          IF( WANTQ ) CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, LI( 1, 1 ), LI( 2, 1 ) );

         // Exit with INFO = 0 if swap was successfully performed.

         RETURN

      } else {

         // CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                 // and 2-by-2 blocks.

         // Solve the generalized Sylvester equation
                  // S11 * R - L * S22 = SCALE * S12
                  // T11 * R - L * T22 = SCALE * T12
         // for R and L. Solutions in LI and IR.

         dlacpy('Full', N1, N2, T( 1, N1+1 ), LDST, LI, LDST );
         dlacpy('Full', N1, N2, S( 1, N1+1 ), LDST, IR( N2+1, N1+1 ), LDST );
         dtgsy2('N', 0, N1, N2, S, LDST, S( N1+1, N1+1 ), LDST, IR( N2+1, N1+1 ), LDST, T, LDST, T( N1+1, N1+1 ), LDST, LI, LDST, SCALE, DSUM, DSCALE, IWORK, IDUM, LINFO );
         if (LINFO != 0) GO TO 70;

         // Compute orthogonal matrix QL:

                     // QL**T * LI = [ TL ]
                                  // [ 0  ]
         // where
                     // LI =  [      -L              ]
                           // [ SCALE * identity(N2) ]

         for (I = 1; I <= N2; I++) { // 10
            dscal(N1, -ONE, LI( 1, I ), 1 );
            LI( N1+I, I ) = SCALE
         } // 10
         dgeqr2(M, N2, LI, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         dorg2r(M, M, N2, LI, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;

         // Compute orthogonal matrix RQ:

                     // IR * RQ**T =   [ 0  TR],

          // where IR = [ SCALE * identity(N1), R ]

         for (I = 1; I <= N1; I++) { // 20
            IR( N2+I, I ) = SCALE
         } // 20
         dgerq2(N1, M, IR( N2+1, 1 ), LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         dorgr2(M, M, N1, IR, LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;

         // Perform the swapping tentatively:

         dgemm('T', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
         dgemm('N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, S, LDST );
         dgemm('T', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
         dgemm('N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, T, LDST );
         dlacpy('F', M, M, S, LDST, SCPY, LDST );
         dlacpy('F', M, M, T, LDST, TCPY, LDST );
         dlacpy('F', M, M, IR, LDST, IRCOP, LDST );
         dlacpy('F', M, M, LI, LDST, LICOP, LDST );

         // Triangularize the B-part by an RQ factorization.
         // Apply transformation (from left) to A-part, giving S.

         dgerq2(M, M, T, LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         CALL DORMR2( 'R', 'T', M, M, M, T, LDST, TAUR, S, LDST, WORK, LINFO )          IF( LINFO != 0 ) GO TO 70;
         CALL DORMR2( 'L', 'N', M, M, M, T, LDST, TAUR, IR, LDST, WORK, LINFO )          IF( LINFO != 0 ) GO TO 70;

         // Compute F-norm(S21) in BRQA21. (T21 is 0.)

         DSCALE = ZERO
         DSUM = ONE
         for (I = 1; I <= N2; I++) { // 30
            dlassq(N1, S( N2+1, I ), 1, DSCALE, DSUM );
         } // 30
         BRQA21 = DSCALE*SQRT( DSUM )

         // Triangularize the B-part by a QR factorization.
         // Apply transformation (from right) to A-part, giving S.

         dgeqr2(M, M, TCPY, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         dorm2r('L', 'T', M, M, M, TCPY, LDST, TAUL, SCPY, LDST, WORK, INFO );
         CALL DORM2R( 'R', 'N', M, M, M, TCPY, LDST, TAUL, LICOP, LDST, WORK, INFO )          IF( LINFO != 0 ) GO TO 70;

         // Compute F-norm(S21) in BQRA21. (T21 is 0.)

         DSCALE = ZERO
         DSUM = ONE
         for (I = 1; I <= N2; I++) { // 40
            dlassq(N1, SCPY( N2+1, I ), 1, DSCALE, DSUM );
         } // 40
         BQRA21 = DSCALE*SQRT( DSUM )

         // Decide which method to use.
           // Weak stability test:
              // F-norm(S21) <= O(EPS * F-norm((S)))

         if ( BQRA21 <= BRQA21 && BQRA21 <= THRESHA ) {
            dlacpy('F', M, M, SCPY, LDST, S, LDST );
            dlacpy('F', M, M, TCPY, LDST, T, LDST );
            dlacpy('F', M, M, IRCOP, LDST, IR, LDST );
            dlacpy('F', M, M, LICOP, LDST, LI, LDST );
         } else if ( BRQA21 >= THRESHA ) {
            GO TO 70
         }

         // Set lower triangle of B-part to zero

         dlaset('Lower', M-1, M-1, ZERO, ZERO, T(2,1), LDST );

         if ( WANDS ) {

            // Strong stability test:
                // F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
                // and
                // F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))

            dlacpy('Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), M );
            dgemm('N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
            dgemm('N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            dlassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SA = DSCALE*SQRT( DSUM )

            dlacpy('Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), M );
            dgemm('N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
            dgemm('N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            dlassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SB = DSCALE*SQRT( DSUM )
            STRONG = SA <= THRESHA && SB <= THRESHB
            if (.NOT.STRONG) GO TO 70;

         }

         // If the swap is accepted ("weakly" and "strongly"), apply the
         // transformations and set N1-by-N2 (2,1)-block to zero.

         dlaset('Full', N1, N2, ZERO, ZERO, S(N2+1,1), LDST );

         // copy back M-by-M diagonal block starting at index J1 of (A, B)

         dlacpy('F', M, M, S, LDST, A( J1, J1 ), LDA );
         dlacpy('F', M, M, T, LDST, B( J1, J1 ), LDB );
         dlaset('Full', LDST, LDST, ZERO, ZERO, T, LDST );

         // Standardize existing 2-by-2 blocks.

         dlaset('Full', M, M, ZERO, ZERO, WORK, M );
         WORK( 1 ) = ONE
         T( 1, 1 ) = ONE
         IDUM = LWORK - M*M - 2
         if ( N2 > 1 ) {
            dlagv2(A( J1, J1 ), LDA, B( J1, J1 ), LDB, AR, AI, BE, WORK( 1 ), WORK( 2 ), T( 1, 1 ), T( 2, 1 ) );
            WORK( M+1 ) = -WORK( 2 )
            WORK( M+2 ) = WORK( 1 )
            T( N2, N2 ) = T( 1, 1 )
            T( 1, 2 ) = -T( 2, 1 )
         }
         WORK( M*M ) = ONE
         T( M, M ) = ONE

         if ( N1 > 1 ) {
            dlagv2(A( J1+N2, J1+N2 ), LDA, B( J1+N2, J1+N2 ), LDB, TAUR, TAUL, WORK( M*M+1 ), WORK( N2*M+N2+1 ), WORK( N2*M+N2+2 ), T( N2+1, N2+1 ), T( M, M-1 ) );
            WORK( M*M ) = WORK( N2*M+N2+1 )
            WORK( M*M-1 ) = -WORK( N2*M+N2+2 )
            T( M, M ) = T( N2+1, N2+1 )
            T( M-1, M ) = -T( M, M-1 )
         }
         dgemm('T', 'N', N2, N1, N2, ONE, WORK, M, A( J1, J1+N2 ), LDA, ZERO, WORK( M*M+1 ), N2 );
         dlacpy('Full', N2, N1, WORK( M*M+1 ), N2, A( J1, J1+N2 ), LDA );
         dgemm('T', 'N', N2, N1, N2, ONE, WORK, M, B( J1, J1+N2 ), LDB, ZERO, WORK( M*M+1 ), N2 );
         dlacpy('Full', N2, N1, WORK( M*M+1 ), N2, B( J1, J1+N2 ), LDB );
         dgemm('N', 'N', M, M, M, ONE, LI, LDST, WORK, M, ZERO, WORK( M*M+1 ), M );
         dlacpy('Full', M, M, WORK( M*M+1 ), M, LI, LDST );
         dgemm('N', 'N', N2, N1, N1, ONE, A( J1, J1+N2 ), LDA, T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 );
         dlacpy('Full', N2, N1, WORK, N2, A( J1, J1+N2 ), LDA );
         dgemm('N', 'N', N2, N1, N1, ONE, B( J1, J1+N2 ), LDB, T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 );
         dlacpy('Full', N2, N1, WORK, N2, B( J1, J1+N2 ), LDB );
         dgemm('T', 'N', M, M, M, ONE, IR, LDST, T, LDST, ZERO, WORK, M );
         dlacpy('Full', M, M, WORK, M, IR, LDST );

         // Accumulate transformations into Q and Z if requested.

         if ( WANTQ ) {
            dgemm('N', 'N', N, M, M, ONE, Q( 1, J1 ), LDQ, LI, LDST, ZERO, WORK, N );
            dlacpy('Full', N, M, WORK, N, Q( 1, J1 ), LDQ );

         }

         if ( WANTZ ) {
            dgemm('N', 'N', N, M, M, ONE, Z( 1, J1 ), LDZ, IR, LDST, ZERO, WORK, N );
            dlacpy('Full', N, M, WORK, N, Z( 1, J1 ), LDZ );

         }

         // Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
                 // (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).

         I = J1 + M
         if ( I <= N ) {
            dgemm('T', 'N', M, N-I+1, M, ONE, LI, LDST, A( J1, I ), LDA, ZERO, WORK, M );
            dlacpy('Full', M, N-I+1, WORK, M, A( J1, I ), LDA );
            dgemm('T', 'N', M, N-I+1, M, ONE, LI, LDST, B( J1, I ), LDB, ZERO, WORK, M );
            dlacpy('Full', M, N-I+1, WORK, M, B( J1, I ), LDB );
         }
         I = J1 - 1
         if ( I > 0 ) {
            dgemm('N', 'N', I, M, M, ONE, A( 1, J1 ), LDA, IR, LDST, ZERO, WORK, I );
            dlacpy('Full', I, M, WORK, I, A( 1, J1 ), LDA );
            dgemm('N', 'N', I, M, M, ONE, B( 1, J1 ), LDB, IR, LDST, ZERO, WORK, I );
            dlacpy('Full', I, M, WORK, I, B( 1, J1 ), LDB );
         }

         // Exit with INFO = 0 if swap was successfully performed.

         RETURN

      }

      // Exit with INFO = 1 if swap was rejected.

      } // 70

      INFO = 1
      RETURN

      // End of DTGEX2

      }
