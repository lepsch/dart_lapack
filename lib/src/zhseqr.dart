      void zhseqr(final int JOB, final int COMPZ, final int N, final int ILO, final int IHI, final Matrix<double> H, final int LDH, final int W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, INFO, LDH, LDZ, LWORK, N;
      String             COMPZ, JOB;
      Complex         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..


      // ==== Matrices of order NTINY or smaller must be processed by
      // .    ZLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== NL allocates some local workspace to help small matrices
      // .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
      // .    required and NL <= NMIN = ilaenv(ISPEC=12,...) is recom-
      // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
      // .    allows up to six simultaneous shifts and a 16-by-16
      // .    deflation window.  ====
      int                NL;
      const              NL = 49 ;
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO;
      const              RZERO = 0.0 ;
      Complex         HL( NL, NL ), WORKL( NL );
      int                KBOT, NMIN;
      bool               INITZ, LQUERY, WANTT, WANTZ;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- bool               lsame;
      // EXTERNAL ILAENV, lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAHQR, ZLAQR0, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN

      // ==== Decode and check the input parameters. ====

      WANTT = lsame( JOB, 'S' );
      INITZ = lsame( COMPZ, 'I' );
      WANTZ = INITZ || lsame( COMPZ, 'V' );
      WORK[1] = DCMPLX( (max( 1, N )).toDouble(), RZERO );
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
         INFO = -10;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO != 0 ) {

         // ==== Quick return in case of invalid argument. ====

         xerbla('ZHSEQR', -INFO );
         return;

      } else if ( N == 0 ) {

         // ==== Quick return in case N = 0; nothing to do. ====

         return;

      } else if ( LQUERY ) {

         // ==== Quick return in case of a workspace query ====

         zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====
         WORK[1] = DCMPLX( max( (WORK( 1 )).toDouble(), (max( 1, N )).toDouble() ), RZERO );
         return;

      } else {

         // ==== copy eigenvalues isolated by ZGEBAL ====

         if (ILO > 1) zcopy( ILO-1, H, LDH+1, W, 1 );
         IF( IHI < N ) zcopy( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 );

         // ==== Initialize Z, if requested ====

         if (INITZ) zlaset( 'A', N, N, ZERO, ONE, Z, LDZ );

         // ==== Quick return if possible ====

         if ( ILO == IHI ) {
            W[ILO] = H( ILO, ILO );
            return;
         }

         // ==== ZLAHQR/ZLAQR0 crossover point ====

         NMIN = ilaenv( 12, 'ZHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK );
         NMIN = max( NTINY, NMIN );

         // ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====

         if ( N > NMIN ) {
            zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );
         } else {

            // ==== Small matrix ====

            zlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO );

            if ( INFO > 0 ) {

               // ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
               // .    when ZLAHQR fails. ====

               KBOT = INFO;

               if ( N >= NL ) {

                  // ==== Larger matrices have enough subdiagonal scratch
                  // .    space to call ZLAQR0 directly. ====

                  zlaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK, INFO );

               } else {

                  // ==== Tiny matrices don't have enough subdiagonal
                  // .    scratch space to benefit from ZLAQR0.  Hence,
                  // .    tiny matrices must be copied into a larger
                  // .    array before calling ZLAQR0. ====

                  zlacpy('A', N, N, H, LDH, HL, NL );
                  HL[N+1][N] = ZERO;
                  zlaset('A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), NL );
                  zlaqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ, WORKL, NL, INFO )                   IF( WANTT || INFO != 0 ) CALL ZLACPY( 'A', N, N, HL, NL, H, LDH );
               }
            }
         }

         // ==== Clear out the trash, if necessary. ====

         if( ( WANTT || INFO != 0 ) && N > 2 ) zlaset( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH );

         // ==== Ensure reported workspace size is backward-compatible with
         // .    previous LAPACK versions. ====

         WORK[1] = DCMPLX( max( (max( 1, N )).toDouble(), (WORK( 1 )).toDouble() ), RZERO );
      }

      // ==== End of ZHSEQR ====

      }
