      void shseqr(JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, LDZ, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      double               H( LDH, * ), WI( * ), WORK( * ), WR( * ), Z( LDZ, * );
      // ..


      // ==== Matrices of order NTINY or smaller must be processed by
      // .    SLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare SLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ilaenv(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               HL( NL, NL ), WORKL( NL );
      int                I, KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLAHQR, SLAQR0, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      // ==== Decode and check the input parameters. ====

      WANTT = lsame( JOB, 'S' );
      INITZ = lsame( COMPZ, 'I' );
      WANTZ = INITZ || lsame( COMPZ, 'V' );
      WORK[1] = SROUNDUP_LWORK( max( 1, N ) );
      LQUERY = LWORK == -1;

      INFO = 0;
      if ( !lsame( JOB, 'E' ) && !WANTT ) {
         INFO = -1;
      } else if ( !lsame( COMPZ, 'N' ) && !WANTZ ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -4;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -5;
      } else if ( LDH < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < max( 1, N ) ) ) {
         INFO = -11;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -13;
      }

      if ( INFO != 0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('SHSEQR', -INFO );
         return;

      } else if ( N == 0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         return;

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         slaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK[1] = max( double( max( 1, N ) ), WORK( 1 ) );
         return;

      } else {

         // ==== copy eigenvalues isolated by SGEBAL ====

         for (I = 1; I <= ILO - 1; I++) { // 10
            WR[I] = H( I, I );
            WI[I] = ZERO;
         } // 10
         for (I = IHI + 1; I <= N; I++) { // 20
            WR[I] = H( I, I );
            WI[I] = ZERO;
         } // 20

         // ==== Initialize Z, if requested ====

         if (INITZ) slaset( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            WR[ILO] = H( ILO, ILO );
            WI[ILO] = ZERO;
            return;
         }

         // ==== SLAHQR/SLAQR0 crossover point ====

         NMIN = ilaenv( 12, 'SHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK );
         NMIN = max( NTINY, NMIN );

         // ==== SLAQR0 for big matrices; SLAHQR for small ones ====

         if ( N > NMIN ) {
            slaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            slahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO );

            if ( INFO > 0 ) {

               // ==== A rare SLAHQR failure!  SLAQR0 sometimes succeeds
               // .    when SLAHQR fails. ====

               KBOT = INFO;

               if ( N >= NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call SLAQR0 directly. ====

                  slaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from SLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling SLAQR0. ====

                  slacpy('A', N, N, H, LDH, HL, NL );
                  HL[N+1][N] = ZERO;
                  slaset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  slaqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT || INFO != 0 ) CALL SLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         if( ( WANTT || INFO != 0 ) && N > 2 ) slaset( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH );

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK[1] = max( double( max( 1, N ) ), WORK( 1 ) );
      }

      // ==== End of SHSEQR ====

      }
