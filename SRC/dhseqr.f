      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      // ..
      // .. Array Arguments ..
      double             H( LDH, * ), WI( * ), WORK( * ), WR( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    DLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare DLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      double             ZERO, ONE;
      const              ZERO = 0.0d0, ONE = 1.0d0 ;
      // ..
      // .. Local Arrays ..
      double             HL( NL, NL ), WORKL( NL );
      // ..
      // .. Local Scalars ..
      int                I, KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      int                ILAENV;
      bool               LSAME;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DLAHQR, DLAQR0, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // ==== Decode and check the input parameters. ====

      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ || LSAME( COMPZ, 'V' )
      WORK( 1 ) = DBLE( MAX( 1, N ) )
      LQUERY = LWORK == -1

      INFO = 0
      if ( .NOT.LSAME( JOB, 'E' ) && .NOT.WANTT ) {
         INFO = -1
      } else if ( .NOT.LSAME( COMPZ, 'N' ) && .NOT.WANTZ ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( ILO < 1 || ILO.GT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( IHI < MIN( ILO, N ) || IHI.GT.N ) {
         INFO = -5
      } else if ( LDH < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDZ < 1 || ( WANTZ && LDZ < MAX( 1, N ) ) ) {
         INFO = -11
      } else if ( LWORK < MAX( 1, N ) && .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO != 0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('DHSEQR', -INFO );
         RETURN

      } else if ( N == 0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         RETURN

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         dlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
         RETURN

      } else {

         // ==== copy eigenvalues isolated by DGEBAL ====

         for (I = 1; I <= ILO - 1; I++) { // 10
            WR( I ) = H( I, I )
            WI( I ) = ZERO
         } // 10
         for (I = IHI + 1; I <= N; I++) { // 20
            WR( I ) = H( I, I )
            WI( I ) = ZERO
         } // 20

         // ==== Initialize Z, if requested ====

         if (INITZ) CALL DLASET( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         }

         // ==== DLAHQR/DLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'DHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== DLAQR0 for big matrices; DLAHQR for small ones ====

         if ( N.GT.NMIN ) {
            dlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            dlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO );

            if ( INFO.GT.0 ) {

               // ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
               // .    when DLAHQR fails. ====

               KBOT = INFO

               if ( N.GE.NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call DLAQR0 directly. ====

                  dlaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from DLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling DLAQR0. ====

                  dlacpy('A', N, N, H, LDH, HL, NL );
                  HL( N+1, N ) = ZERO
                  dlaset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  dlaqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT || INFO != 0 ) CALL DLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         IF( ( WANTT || INFO != 0 ) && N.GT.2 ) CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
      }

      // ==== End of DHSEQR ====

      }
