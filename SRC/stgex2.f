      SUBROUTINE STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, N1, N2, WORK, LWORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================
*  Replaced various illegal calls to SCOPY by calls to SLASET, or by DO
*  loops. Sven Hammarling, 1/5/02.

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               TWENTY
      const              TWENTY = 2.0E+01 ;
      int                LDST;
      const              LDST = 4 ;
      bool               WANDS;
      const              WANDS = true ;
      // ..
      // .. Local Scalars ..
      bool               STRONG, WEAK;
      int                I, IDUM, LINFO, M;
      REAL               BQRA21, BRQA21, DDUM, DNORMA, DNORMB, DSCALE, DSUM, EPS, F, G, SA, SB, SCALE, SMLNUM, THRESHA, THRESHB
      // ..
      // .. Local Arrays ..
      int                IWORK( LDST + 2 );
      REAL               AI( 2 ), AR( 2 ), BE( 2 ), IR( LDST, LDST ), IRCOP( LDST, LDST ), LI( LDST, LDST ), LICOP( LDST, LDST ), S( LDST, LDST ), SCPY( LDST, LDST ), T( LDST, LDST ), TAUL( LDST ), TAUR( LDST ), TCPY( LDST, LDST )
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGEQR2, SGERQ2, SLACPY, SLAGV2, SLARTG, SLASET, SLASSQ, SORG2R, SORGR2, SORM2R, SORMR2, SROT, SSCAL, STGSY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if (N <= 1 || N1 <= 0 || N2 <= 0) RETURN       IF( N1 > N || ( J1+N1 ) > N ) RETURN;
      M = N1 + N2
      if ( LWORK < MAX( N*M, M*M*2 ) ) {
         INFO = -16
         WORK( 1 ) = MAX( N*M, M*M*2 )
         RETURN
      }

      WEAK = false;
      STRONG = false;

      // Make a local copy of selected block

      slaset('Full', LDST, LDST, ZERO, ZERO, LI, LDST );
      slaset('Full', LDST, LDST, ZERO, ZERO, IR, LDST );
      slacpy('Full', M, M, A( J1, J1 ), LDA, S, LDST );
      slacpy('Full', M, M, B( J1, J1 ), LDB, T, LDST );

      // Compute threshold for testing acceptance of swapping.

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      DSCALE = ZERO
      DSUM = ONE
      slacpy('Full', M, M, S, LDST, WORK, M );
      slassq(M*M, WORK, 1, DSCALE, DSUM );
      DNORMA = DSCALE*SQRT( DSUM )
      DSCALE = ZERO
      DSUM = ONE
      slacpy('Full', M, M, T, LDST, WORK, M );
      slassq(M*M, WORK, 1, DSCALE, DSUM );
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
         slartg(F, G, IR( 1, 2 ), IR( 1, 1 ), DDUM );
         IR( 2, 1 ) = -IR( 1, 2 )
         IR( 2, 2 ) = IR( 1, 1 )
         srot(2, S( 1, 1 ), 1, S( 1, 2 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         srot(2, T( 1, 1 ), 1, T( 1, 2 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         if ( SA >= SB ) {
            slartg(S( 1, 1 ), S( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), DDUM );
         } else {
            slartg(T( 1, 1 ), T( 2, 1 ), LI( 1, 1 ), LI( 2, 1 ), DDUM );
         }
         srot(2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, LI( 1, 1 ), LI( 2, 1 ) );
         srot(2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, LI( 1, 1 ), LI( 2, 1 ) );
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

            slacpy('Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), M );
            sgemm('N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
            sgemm('N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            slassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SA = DSCALE*SQRT( DSUM )

            slacpy('Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), M );
            sgemm('N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
            sgemm('N', 'T', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            slassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SB = DSCALE*SQRT( DSUM )
            STRONG = SA <= THRESHA && SB <= THRESHB
            if (.NOT.STRONG) GO TO 70;
         }

         // Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
                // (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).

         srot(J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         srot(J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) );
         srot(N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, LI( 1, 1 ), LI( 2, 1 ) );
         srot(N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, LI( 1, 1 ), LI( 2, 1 ) );

         // Set  N1-by-N2 (2,1) - blocks to ZERO.

         A( J1+1, J1 ) = ZERO
         B( J1+1, J1 ) = ZERO

         // Accumulate transformations into Q and Z if requested.

         if (WANTZ) CALL SROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, IR( 1, 1 ), IR( 2, 1 ) )          IF( WANTQ ) CALL SROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, LI( 1, 1 ), LI( 2, 1 ) );

         // Exit with INFO = 0 if swap was successfully performed.

         RETURN

      } else {

         // CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
                 // and 2-by-2 blocks.

         // Solve the generalized Sylvester equation
                  // S11 * R - L * S22 = SCALE * S12
                  // T11 * R - L * T22 = SCALE * T12
         // for R and L. Solutions in LI and IR.

         slacpy('Full', N1, N2, T( 1, N1+1 ), LDST, LI, LDST );
         slacpy('Full', N1, N2, S( 1, N1+1 ), LDST, IR( N2+1, N1+1 ), LDST );
         stgsy2('N', 0, N1, N2, S, LDST, S( N1+1, N1+1 ), LDST, IR( N2+1, N1+1 ), LDST, T, LDST, T( N1+1, N1+1 ), LDST, LI, LDST, SCALE, DSUM, DSCALE, IWORK, IDUM, LINFO );
         if (LINFO != 0) GO TO 70;

         // Compute orthogonal matrix QL:

                     // QL**T * LI = [ TL ]
                                  // [ 0  ]
         // where
                     // LI =  [      -L              ]
                           // [ SCALE * identity(N2) ]

         for (I = 1; I <= N2; I++) { // 10
            sscal(N1, -ONE, LI( 1, I ), 1 );
            LI( N1+I, I ) = SCALE
         } // 10
         sgeqr2(M, N2, LI, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         sorg2r(M, M, N2, LI, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;

         // Compute orthogonal matrix RQ:

                     // IR * RQ**T =   [ 0  TR],

          // where IR = [ SCALE * identity(N1), R ]

         for (I = 1; I <= N1; I++) { // 20
            IR( N2+I, I ) = SCALE
         } // 20
         sgerq2(N1, M, IR( N2+1, 1 ), LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         sorgr2(M, M, N1, IR, LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;

         // Perform the swapping tentatively:

         sgemm('T', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
         sgemm('N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, S, LDST );
         sgemm('T', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
         sgemm('N', 'T', M, M, M, ONE, WORK, M, IR, LDST, ZERO, T, LDST );
         slacpy('F', M, M, S, LDST, SCPY, LDST );
         slacpy('F', M, M, T, LDST, TCPY, LDST );
         slacpy('F', M, M, IR, LDST, IRCOP, LDST );
         slacpy('F', M, M, LI, LDST, LICOP, LDST );

         // Triangularize the B-part by an RQ factorization.
         // Apply transformation (from left) to A-part, giving S.

         sgerq2(M, M, T, LDST, TAUR, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         CALL SORMR2( 'R', 'T', M, M, M, T, LDST, TAUR, S, LDST, WORK, LINFO )          IF( LINFO != 0 ) GO TO 70;
         CALL SORMR2( 'L', 'N', M, M, M, T, LDST, TAUR, IR, LDST, WORK, LINFO )          IF( LINFO != 0 ) GO TO 70;

         // Compute F-norm(S21) in BRQA21. (T21 is 0.)

         DSCALE = ZERO
         DSUM = ONE
         for (I = 1; I <= N2; I++) { // 30
            slassq(N1, S( N2+1, I ), 1, DSCALE, DSUM );
         } // 30
         BRQA21 = DSCALE*SQRT( DSUM )

         // Triangularize the B-part by a QR factorization.
         // Apply transformation (from right) to A-part, giving S.

         sgeqr2(M, M, TCPY, LDST, TAUL, WORK, LINFO );
         if (LINFO != 0) GO TO 70;
         sorm2r('L', 'T', M, M, M, TCPY, LDST, TAUL, SCPY, LDST, WORK, INFO );
         CALL SORM2R( 'R', 'N', M, M, M, TCPY, LDST, TAUL, LICOP, LDST, WORK, INFO )          IF( LINFO != 0 ) GO TO 70;

         // Compute F-norm(S21) in BQRA21. (T21 is 0.)

         DSCALE = ZERO
         DSUM = ONE
         for (I = 1; I <= N2; I++) { // 40
            slassq(N1, SCPY( N2+1, I ), 1, DSCALE, DSUM );
         } // 40
         BQRA21 = DSCALE*SQRT( DSUM )

         // Decide which method to use.
           // Weak stability test:
              // F-norm(S21) <= O(EPS * F-norm((S)))

         if ( BQRA21 <= BRQA21 && BQRA21 <= THRESHA ) {
            slacpy('F', M, M, SCPY, LDST, S, LDST );
            slacpy('F', M, M, TCPY, LDST, T, LDST );
            slacpy('F', M, M, IRCOP, LDST, IR, LDST );
            slacpy('F', M, M, LICOP, LDST, LI, LDST );
         } else if ( BRQA21 >= THRESHA ) {
            GO TO 70
         }

         // Set lower triangle of B-part to zero

         slaset('Lower', M-1, M-1, ZERO, ZERO, T(2,1), LDST );

         if ( WANDS ) {

            // Strong stability test:
                // F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
                // and
                // F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))

            slacpy('Full', M, M, A( J1, J1 ), LDA, WORK( M*M+1 ), M );
            sgemm('N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK, M );
            sgemm('N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            slassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SA = DSCALE*SQRT( DSUM )

            slacpy('Full', M, M, B( J1, J1 ), LDB, WORK( M*M+1 ), M );
            sgemm('N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK, M );
            sgemm('N', 'N', M, M, M, -ONE, WORK, M, IR, LDST, ONE, WORK( M*M+1 ), M );
            DSCALE = ZERO
            DSUM = ONE
            slassq(M*M, WORK( M*M+1 ), 1, DSCALE, DSUM );
            SB = DSCALE*SQRT( DSUM )
            STRONG = SA <= THRESHA && SB <= THRESHB
            if (.NOT.STRONG) GO TO 70;

         }

         // If the swap is accepted ("weakly" and "strongly"), apply the
         // transformations and set N1-by-N2 (2,1)-block to zero.

         slaset('Full', N1, N2, ZERO, ZERO, S(N2+1,1), LDST );

         // copy back M-by-M diagonal block starting at index J1 of (A, B)

         slacpy('F', M, M, S, LDST, A( J1, J1 ), LDA );
         slacpy('F', M, M, T, LDST, B( J1, J1 ), LDB );
         slaset('Full', LDST, LDST, ZERO, ZERO, T, LDST );

         // Standardize existing 2-by-2 blocks.

         slaset('Full', M, M, ZERO, ZERO, WORK, M );
         WORK( 1 ) = ONE
         T( 1, 1 ) = ONE
         IDUM = LWORK - M*M - 2
         if ( N2 > 1 ) {
            slagv2(A( J1, J1 ), LDA, B( J1, J1 ), LDB, AR, AI, BE, WORK( 1 ), WORK( 2 ), T( 1, 1 ), T( 2, 1 ) );
            WORK( M+1 ) = -WORK( 2 )
            WORK( M+2 ) = WORK( 1 )
            T( N2, N2 ) = T( 1, 1 )
            T( 1, 2 ) = -T( 2, 1 )
         }
         WORK( M*M ) = ONE
         T( M, M ) = ONE

         if ( N1 > 1 ) {
            slagv2(A( J1+N2, J1+N2 ), LDA, B( J1+N2, J1+N2 ), LDB, TAUR, TAUL, WORK( M*M+1 ), WORK( N2*M+N2+1 ), WORK( N2*M+N2+2 ), T( N2+1, N2+1 ), T( M, M-1 ) );
            WORK( M*M ) = WORK( N2*M+N2+1 )
            WORK( M*M-1 ) = -WORK( N2*M+N2+2 )
            T( M, M ) = T( N2+1, N2+1 )
            T( M-1, M ) = -T( M, M-1 )
         }
         sgemm('T', 'N', N2, N1, N2, ONE, WORK, M, A( J1, J1+N2 ), LDA, ZERO, WORK( M*M+1 ), N2 );
         slacpy('Full', N2, N1, WORK( M*M+1 ), N2, A( J1, J1+N2 ), LDA );
         sgemm('T', 'N', N2, N1, N2, ONE, WORK, M, B( J1, J1+N2 ), LDB, ZERO, WORK( M*M+1 ), N2 );
         slacpy('Full', N2, N1, WORK( M*M+1 ), N2, B( J1, J1+N2 ), LDB );
         sgemm('N', 'N', M, M, M, ONE, LI, LDST, WORK, M, ZERO, WORK( M*M+1 ), M );
         slacpy('Full', M, M, WORK( M*M+1 ), M, LI, LDST );
         sgemm('N', 'N', N2, N1, N1, ONE, A( J1, J1+N2 ), LDA, T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 );
         slacpy('Full', N2, N1, WORK, N2, A( J1, J1+N2 ), LDA );
         sgemm('N', 'N', N2, N1, N1, ONE, B( J1, J1+N2 ), LDB, T( N2+1, N2+1 ), LDST, ZERO, WORK, N2 );
         slacpy('Full', N2, N1, WORK, N2, B( J1, J1+N2 ), LDB );
         sgemm('T', 'N', M, M, M, ONE, IR, LDST, T, LDST, ZERO, WORK, M );
         slacpy('Full', M, M, WORK, M, IR, LDST );

         // Accumulate transformations into Q and Z if requested.

         if ( WANTQ ) {
            sgemm('N', 'N', N, M, M, ONE, Q( 1, J1 ), LDQ, LI, LDST, ZERO, WORK, N );
            slacpy('Full', N, M, WORK, N, Q( 1, J1 ), LDQ );

         }

         if ( WANTZ ) {
            sgemm('N', 'N', N, M, M, ONE, Z( 1, J1 ), LDZ, IR, LDST, ZERO, WORK, N );
            slacpy('Full', N, M, WORK, N, Z( 1, J1 ), LDZ );

         }

         // Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
                 // (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).

         I = J1 + M
         if ( I <= N ) {
            sgemm('T', 'N', M, N-I+1, M, ONE, LI, LDST, A( J1, I ), LDA, ZERO, WORK, M );
            slacpy('Full', M, N-I+1, WORK, M, A( J1, I ), LDA );
            sgemm('T', 'N', M, N-I+1, M, ONE, LI, LDST, B( J1, I ), LDB, ZERO, WORK, M );
            slacpy('Full', M, N-I+1, WORK, M, B( J1, I ), LDB );
         }
         I = J1 - 1
         if ( I > 0 ) {
            sgemm('N', 'N', I, M, M, ONE, A( 1, J1 ), LDA, IR, LDST, ZERO, WORK, I );
            slacpy('Full', I, M, WORK, I, A( 1, J1 ), LDA );
            sgemm('N', 'N', I, M, M, ONE, B( 1, J1 ), LDB, IR, LDST, ZERO, WORK, I );
            slacpy('Full', I, M, WORK, I, B( 1, J1 ), LDB );
         }

         // Exit with INFO = 0 if swap was successfully performed.

         RETURN

      }

      // Exit with INFO = 1 if swap was rejected.

      } // 70

      INFO = 1
      RETURN

      // End of STGEX2

      }
