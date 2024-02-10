      void zlaqz3(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, NBLOCK_DESIRED, ALPHA, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, QC, LDQC, ZC, LDZC, WORK, LWORK, Box<int> INFO ) {
      // Function arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC;
       Complex, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), ALPHA( * ), BETA( * );

      int    , INTENT( OUT ) :: INFO;

      // Parameters
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      int     :: I, J, NS, ISTARTM, ISTOPM, SHEIGHT, SWIDTH, K, NP, ISTARTB, ISTOPB, ISHIFT, NBLOCK, NPOS;
      double           :: SAFMIN, SAFMAX, C, SCALE;
      Complex :: TEMP, TEMP2, TEMP3, S;

      // External Functions
      // EXTERNAL :: XERBLA, ZLASET, ZLARTG, ZROT, ZLAQZ1, ZGEMM, ZLACPY
      double          , EXTERNAL :: DLAMCH;

      INFO = 0;
      if ( NBLOCK_DESIRED < NSHIFTS+1 ) {
         INFO = -8;
      }
      if ( LWORK == -1 ) {
         // workspace query, quick return;
         WORK[1] = N*NBLOCK_DESIRED;
         return;
      } else if ( LWORK < N*NBLOCK_DESIRED ) {
         INFO = -25;
      }

      if ( INFO != 0 ) {
         xerbla('ZLAQZ3', -INFO );
         return;
      }


      // Executable statements


      // Get machine constants
      SAFMIN = dlamch( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;

      if ( ILO >= IHI ) {
         return;
      }

      if ( ILSCHUR ) {
         ISTARTM = 1;
         ISTOPM = N;
      } else {
         ISTARTM = ILO;
         ISTOPM = IHI;
      }

      NS = NSHIFTS;
      NPOS = max( NBLOCK_DESIRED-NS, 1 );


      // The following block introduces the shifts and chases
      // them down one by one just enough to make space for
      // the other shifts. The near-the-diagonal block is
      // of size (ns+1) x ns.

      zlaset('FULL', NS+1, NS+1, CZERO, CONE, QC, LDQC );
      zlaset('FULL', NS, NS, CZERO, CONE, ZC, LDZC );

      for (I = 1; I <= NS; I++) {
         // Introduce the shift
         SCALE = sqrt( ( ALPHA( I ) ).abs() ) * sqrt( ( BETA( I ) ).abs() );
         if ( SCALE >= SAFMIN && SCALE <= SAFMAX ) {
            ALPHA[I] = ALPHA( I )/SCALE;
            BETA[I] = BETA( I )/SCALE;
         }

         TEMP2 = BETA( I )*A( ILO, ILO )-ALPHA( I )*B( ILO, ILO );
         TEMP3 = BETA( I )*A( ILO+1, ILO );
          if ( ( TEMP2 ).abs() > SAFMAX || ( TEMP3 ).abs() > SAFMAX ) {
            TEMP2 = CONE;
            TEMP3 = CZERO;
         }

         zlartg(TEMP2, TEMP3, C, S, TEMP );
         zrot(NS, A( ILO, ILO ), LDA, A( ILO+1, ILO ), LDA, C, S );
         zrot(NS, B( ILO, ILO ), LDB, B( ILO+1, ILO ), LDB, C, S );
         zrot(NS+1, QC( 1, 1 ), 1, QC( 1, 2 ), 1, C, DCONJG( S ) );

         // Chase the shift down
         for (J = 1; J <= NS-I; J++) {
             zlaqz1( true , true , J, 1, NS, IHI-ILO+1, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, NS+1, 1, QC, LDQC, NS, 1, ZC, LDZC );

         }

      }

      // Update the rest of the pencil

      // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
      // from the left with Qc(1:ns+1,1:ns+1)'
      SHEIGHT = NS+1;
      SWIDTH = ISTOPM-( ILO+NS )+1;
      if ( SWIDTH > 0 ) {
         zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( ILO, ILO+NS ), LDA, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ILO, ILO+NS ), LDA );
         zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( ILO, ILO+NS ), LDB, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ILO, ILO+NS ), LDB );
      }
      if ( ILQ ) {
         zgemm('N', 'N', N, SHEIGHT, SHEIGHT, CONE, Q( 1, ILO ), LDQ, QC, LDQC, CZERO, WORK, N );
         zlacpy('ALL', N, SHEIGHT, WORK, N, Q( 1, ILO ), LDQ );
      }

      // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
      // from the right with Zc(1:ns,1:ns)
      SHEIGHT = ILO-1-ISTARTM+1;
      SWIDTH = NS;
      if ( SHEIGHT > 0 ) {
         zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, ILO ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, ILO ), LDA );
         zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, ILO ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, ILO ), LDB );
      }
      if ( ILZ ) {
         zgemm('N', 'N', N, SWIDTH, SWIDTH, CONE, Z( 1, ILO ), LDZ, ZC, LDZC, CZERO, WORK, N );
         zlacpy('ALL', N, SWIDTH, WORK, N, Z( 1, ILO ), LDZ );
      }

      // The following block chases the shifts down to the bottom
      // right block. If possible, a shift is moved down npos
      // positions at a time

      K = ILO;
      while (K < IHI-NS) {
         NP = min( IHI-NS-K, NPOS );
         // Size of the near-the-diagonal block
         NBLOCK = NS+NP;
         // istartb points to the first row we will be updating
         ISTARTB = K+1;
         // istopb points to the last column we will be updating
         ISTOPB = K+NBLOCK-1;

         zlaset('FULL', NS+NP, NS+NP, CZERO, CONE, QC, LDQC );
         zlaset('FULL', NS+NP, NS+NP, CZERO, CONE, ZC, LDZC );

         // Near the diagonal shift chase
         for (I = NS-1; I >= 0; I--) {
            for (J = 0; J <= NP-1; J++) {
               // Move down the block with index k+i+j, updating
               // the (ns+np x ns+np) block:
               // (k:k+ns+np,k:k+ns+np-1)
               zlaqz1( true , true , K+I+J, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NBLOCK, K+1, QC, LDQC, NBLOCK, K, ZC, LDZC );
            }
         }

         // Update rest of the pencil

         // Update A(k+1:k+ns+np, k+ns+np:istopm) and
         // B(k+1:k+ns+np, k+ns+np:istopm)
         // from the left with Qc(1:ns+np,1:ns+np)'
         SHEIGHT = NS+NP;
         SWIDTH = ISTOPM-( K+NS+NP )+1;
         if ( SWIDTH > 0 ) {
            zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( K+1, K+NS+NP ), LDA, CZERO, WORK, SHEIGHT );
            zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( K+1, K+NS+NP ), LDA );
            zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( K+1, K+NS+NP ), LDB, CZERO, WORK, SHEIGHT );
            zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( K+1, K+NS+NP ), LDB );
         }
         if ( ILQ ) {
            zgemm('N', 'N', N, NBLOCK, NBLOCK, CONE, Q( 1, K+1 ), LDQ, QC, LDQC, CZERO, WORK, N );
            zlacpy('ALL', N, NBLOCK, WORK, N, Q( 1, K+1 ), LDQ );
         }

         // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
         // from the right with Zc(1:ns+np,1:ns+np)
         SHEIGHT = K-ISTARTM+1;
         SWIDTH = NBLOCK;
         if ( SHEIGHT > 0 ) {
            zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, K ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT );
            zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, K ), LDA );
            zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, K ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT );
            zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, K ), LDB );
         }
         if ( ILZ ) {
            zgemm('N', 'N', N, NBLOCK, NBLOCK, CONE, Z( 1, K ), LDZ, ZC, LDZC, CZERO, WORK, N );
            zlacpy('ALL', N, NBLOCK, WORK, N, Z( 1, K ), LDZ );
         }

         K = K+NP;

      }

      // The following block removes the shifts from the bottom right corner
      // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      zlaset('FULL', NS, NS, CZERO, CONE, QC, LDQC );
      zlaset('FULL', NS+1, NS+1, CZERO, CONE, ZC, LDZC );

      // istartb points to the first row we will be updating
      ISTARTB = IHI-NS+1;
      // istopb points to the last column we will be updating
      ISTOPB = IHI;

      for (I = 1; I <= NS; I++) {
         // Chase the shift down to the bottom right corner
         for (ISHIFT = IHI-I; ISHIFT <= IHI-1; ISHIFT++) {
            zlaqz1( true , true , ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS, IHI-NS+1, QC, LDQC, NS+1, IHI-NS, ZC, LDZC );
         }

      }

      // Update rest of the pencil

      // Update A(ihi-ns+1:ihi, ihi+1:istopm)
      // from the left with Qc(1:ns,1:ns)'
      SHEIGHT = NS;
      SWIDTH = ISTOPM-( IHI+1 )+1;
      if ( SWIDTH > 0 ) {
         zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( IHI-NS+1, IHI+1 ), LDA, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( IHI-NS+1, IHI+1 ), LDA );
         zgemm('C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( IHI-NS+1, IHI+1 ), LDB, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( IHI-NS+1, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         zgemm('N', 'N', N, NS, NS, CONE, Q( 1, IHI-NS+1 ), LDQ, QC, LDQC, CZERO, WORK, N );
         zlacpy('ALL', N, NS, WORK, N, Q( 1, IHI-NS+1 ), LDQ );
      }

      // Update A(istartm:ihi-ns,ihi-ns:ihi)
      // from the right with Zc(1:ns+1,1:ns+1)
      SHEIGHT = IHI-NS-ISTARTM+1;
      SWIDTH = NS+1;
      if ( SHEIGHT > 0 ) {
         zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, IHI-NS ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, IHI-NS ), LDA );
         zgemm('N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, IHI-NS ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT );
         zlacpy('ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, IHI-NS ), LDB );
      }
      if ( ILZ ) {
         zgemm('N', 'N', N, NS+1, NS+1, CONE, Z( 1, IHI-NS ), LDZ, ZC, LDZC, CZERO, WORK, N );
         zlacpy('ALL', N, NS+1, WORK, N, Z( 1, IHI-NS ), LDZ );
      }

      END SUBROUTINE;
