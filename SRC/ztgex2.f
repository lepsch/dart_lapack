      SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                INFO, J1, LDA, LDB, LDQ, LDZ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      double             TWENTY;
      const              TWENTY = 2.0D+1 ;
      int                LDST;
      const              LDST = 2 ;
      bool               WANDS;
      const              WANDS = true ;
      // ..
      // .. Local Scalars ..
      bool               STRONG, WEAK;
      int                I, M;
      double             CQ, CZ, EPS, SA, SB, SCALE, SMLNUM, SUM, THRESHA, THRESHB;
      COMPLEX*16         CDUM, F, G, SQ, SZ
      // ..
      // .. Local Arrays ..
      COMPLEX*16         S( LDST, LDST ), T( LDST, LDST ), WORK( 8 )
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACPY, ZLARTG, ZLASSQ, ZROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if (N.LE.1) RETURN;

      M = LDST
      WEAK = false;
      STRONG = false;

      // Make a local copy of selected block in (A, B)

      zlacpy('Full', M, M, A( J1, J1 ), LDA, S, LDST );
      zlacpy('Full', M, M, B( J1, J1 ), LDB, T, LDST );

      // Compute the threshold for testing the acceptance of swapping.

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SCALE = DBLE( CZERO )
      SUM = DBLE( CONE )
      zlacpy('Full', M, M, S, LDST, WORK, M );
      zlacpy('Full', M, M, T, LDST, WORK( M*M+1 ), M );
      zlassq(M*M, WORK, 1, SCALE, SUM );
      SA = SCALE*SQRT( SUM )
      SCALE = DBLE( CZERO )
      SUM = DBLE( CONE )
      zlassq(M*M, WORK(M*M+1), 1, SCALE, SUM );
      SB = SCALE*SQRT( SUM )

      // THRES has been changed from
         // THRESH = MAX( TEN*EPS*SA, SMLNUM )
      // to
         // THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
      // on 04/01/10.
      // "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
      // Jim Demmel and Guillaume Revy. See forum post 1783.

      THRESHA = MAX( TWENTY*EPS*SA, SMLNUM )
      THRESHB = MAX( TWENTY*EPS*SB, SMLNUM )

      // Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
      // using Givens rotations and perform the swap tentatively.

      F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
      G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
      SA = ABS( S( 2, 2 ) ) * ABS( T( 1, 1 ) )
      SB = ABS( S( 1, 1 ) ) * ABS( T( 2, 2 ) )
      zlartg(G, F, CZ, SZ, CDUM );
      SZ = -SZ
      zrot(2, S( 1, 1 ), 1, S( 1, 2 ), 1, CZ, DCONJG( SZ ) );
      zrot(2, T( 1, 1 ), 1, T( 1, 2 ), 1, CZ, DCONJG( SZ ) );
      if ( SA >= SB ) {
         zlartg(S( 1, 1 ), S( 2, 1 ), CQ, SQ, CDUM );
      } else {
         zlartg(T( 1, 1 ), T( 2, 1 ), CQ, SQ, CDUM );
      }
      zrot(2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, CQ, SQ );
      zrot(2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, CQ, SQ );

      // Weak stability test: |S21| <= O(EPS F-norm((A)))
                           // and  |T21| <= O(EPS F-norm((B)))

      WEAK = ABS( S( 2, 1 ) ).LE.THRESHA && ABS( T( 2, 1 ) ).LE.THRESHB       IF( .NOT.WEAK ) GO TO 20

      if ( WANDS ) {

         // Strong stability test:
            // F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
            // and
            // F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))

         zlacpy('Full', M, M, S, LDST, WORK, M );
         zlacpy('Full', M, M, T, LDST, WORK( M*M+1 ), M );
         zrot(2, WORK, 1, WORK( 3 ), 1, CZ, -DCONJG( SZ ) );
         zrot(2, WORK( 5 ), 1, WORK( 7 ), 1, CZ, -DCONJG( SZ ) );
         zrot(2, WORK, 2, WORK( 2 ), 2, CQ, -SQ );
         zrot(2, WORK( 5 ), 2, WORK( 6 ), 2, CQ, -SQ );
         for (I = 1; I <= 2; I++) { // 10
            WORK( I ) = WORK( I ) - A( J1+I-1, J1 )
            WORK( I+2 ) = WORK( I+2 ) - A( J1+I-1, J1+1 )
            WORK( I+4 ) = WORK( I+4 ) - B( J1+I-1, J1 )
            WORK( I+6 ) = WORK( I+6 ) - B( J1+I-1, J1+1 )
         } // 10
         SCALE = DBLE( CZERO )
         SUM = DBLE( CONE )
         zlassq(M*M, WORK, 1, SCALE, SUM );
         SA = SCALE*SQRT( SUM )
         SCALE = DBLE( CZERO )
         SUM = DBLE( CONE )
         zlassq(M*M, WORK(M*M+1), 1, SCALE, SUM );
         SB = SCALE*SQRT( SUM )
         STRONG = SA.LE.THRESHA && SB.LE.THRESHB
         if (.NOT.STRONG) GO TO 20;
      }

      // If the swap is accepted ("weakly" and "strongly"), apply the
      // equivalence transformations to the original matrix pair (A,B)

      zrot(J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, CZ, DCONJG( SZ ) )       CALL ZROT( J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, CZ, DCONJG( SZ ) );
      zrot(N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, CQ, SQ );
      zrot(N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, CQ, SQ );

      // Set  N1 by N2 (2,1) blocks to 0

      A( J1+1, J1 ) = CZERO
      B( J1+1, J1 ) = CZERO

      // Accumulate transformations into Q and Z if requested.

      if (WANTZ) CALL ZROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, CZ, DCONJG( SZ ) )       IF( WANTQ ) CALL ZROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, CQ, DCONJG( SQ ) );

      // Exit with INFO = 0 if swap was successfully performed.

      RETURN

      // Exit with INFO = 1 if swap was rejected.

      } // 20
      INFO = 1
      RETURN

      // End of ZTGEX2

      }
